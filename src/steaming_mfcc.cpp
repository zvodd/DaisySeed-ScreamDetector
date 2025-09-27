#include "daisy_seed.h"
#include "daisysp.h"
#include <vector>
#include <string>
#include <cmath>

#include "arm_math.h"



using namespace daisy;
using namespace daisysp;

// ====================================================================
// 1. Configuration Constants
// ====================================================================

// Audio settings
const float SAMPLE_RATE = 16000.0f; // Lower sample rate for voice
const int   BLOCK_SIZE  = 64;       // Number of samples per audio callback

// MFCC Frame processing settings
const int   FRAME_LENGTH      = 400;  // Window size in samples (25ms * 16000Hz)
const int   FRAME_STRIDE      = 160;  // Hop size in samples (10ms * 16000Hz)
const int   FFT_SIZE          = 512;  // Must be >= FRAME_LENGTH and a power of 2
const float PRE_EMPHASIS_ALPHA = 0.97f;

// Mel Filterbank settings
const int NUM_MEL_BANDS     = 40;
const int NUM_MFCC_COEFFS   = 13;   // Number of coefficients to keep

// Circular buffer to hold incoming audio
// Should be large enough to prevent overflow
const int CIRCULAR_BUFFER_SIZE = 2048; 

// ====================================================================
// 2. Global Buffers & Variables
// ====================================================================

DaisySeed hw;

// Circular buffer for audio input
float DSY_SDRAM_BSS circular_audio_buffer[CIRCULAR_BUFFER_SIZE];
volatile int write_pos = 0;
int read_pos = 0;

// MFCC processing buffers
float frame[FRAME_LENGTH];
float fft_input[FFT_SIZE];
float fft_output[FFT_SIZE];
float power_spectrum[FFT_SIZE / 2 + 1];
float mel_energies[NUM_MEL_BANDS];
float mfcc_coeffs[NUM_MFCC_COEFFS];

// Pre-calculated data
float hamming_window[FRAME_LENGTH];
std::vector<std::vector<float>> mel_filterbank;
float prev_sample_for_preemphasis = 0.0f; // Store last sample for pre-emphasis continuity

// CMSIS-DSP FFT instance
arm_rfft_fast_instance_f32 fft_instance;

// ====================================================================
// 3. Helper and Initialization Functions
// ====================================================================

// Hz to Mel scale conversion
inline float hz_to_mel(float hz) {
    return 2595.0f * log10f(1.0f + hz / 700.0f);
}

// Mel to Hz scale conversion
inline float mel_to_hz(float mel) {
    return 700.0f * (powf(10.0f, mel / 2595.0f) - 1.0f);
}

/**
 * @brief Basic DCT-II implementation.
 * @note For production, consider replacing the entire MFCC chain with 
 * the highly optimized arm_mfcc_f32 from CMSIS-DSP.
 */
void dct_ii(float* input, float* output, int size, int num_coeffs) {
    for (int k = 0; k < num_coeffs; ++k) {
        float sum = 0.0f;
        for (int n = 0; n < size; ++n) {
            sum += input[n] * cosf(M_PI / size * (n + 0.5f) * k);
        }
        float scale = (k == 0) ? sqrtf(1.0f / size) : sqrtf(2.0f / size);
        output[k] = sum * scale; // Apply ortho-normal scaling
    }
}

void init_hamming_window() {
    for (size_t i = 0; i < FRAME_LENGTH; ++i) {
        hamming_window[i] = 0.54f - 0.46f * cosf(2.0f * M_PI * i / (FRAME_LENGTH - 1));
    }
}

// (init_mel_filterbank function is unchanged from your original code)
void init_mel_filterbank() {
    mel_filterbank.resize(NUM_MEL_BANDS, std::vector<float>(FFT_SIZE / 2 + 1, 0.0f));

    float low_freq_mel = 0;
    float high_freq_mel = hz_to_mel(SAMPLE_RATE / 2.0f);
    std::vector<float> mel_points(NUM_MEL_BANDS + 2);
    for (int i = 0; i < NUM_MEL_BANDS + 2; ++i) {
        mel_points[i] = low_freq_mel + i * (high_freq_mel - low_freq_mel) / (NUM_MEL_BANDS + 1);
    }

    std::vector<float> hz_points(NUM_MEL_BANDS + 2);
    std::vector<size_t> fft_bins(NUM_MEL_BANDS + 2);
    for (int i = 0; i < NUM_MEL_BANDS + 2; ++i) {
        hz_points[i] = mel_to_hz(mel_points[i]);
        fft_bins[i] = static_cast<size_t>(floorf((FFT_SIZE + 1) * hz_points[i] / SAMPLE_RATE));
    }

    for (int m = 1; m < NUM_MEL_BANDS + 1; ++m) {
        size_t f_m_minus = fft_bins[m - 1];
        size_t f_m = fft_bins[m];
        size_t f_m_plus = fft_bins[m + 1];

        for (size_t k = f_m_minus; k < f_m; ++k) {
            mel_filterbank[m - 1][k] = (float)(k - f_m_minus) / (f_m - f_m_minus);
        }
        for (size_t k = f_m; k < f_m_plus; ++k) {
            mel_filterbank[m - 1][k] = (float)(f_m_plus - k) / (f_m_plus - f_m);
        }
    }
}

// We need to declare the external table struct we want to use
extern arm_cfft_instance_f32 arm_cfft_sR_f32_len2048;
extern const float32_t twiddleCoef_rfft_4096[];

/**
 * @brief  Specialized RFFT init function for a fixed size of 2048.
 * @param[in,out] S  points to an arm_rfft_fast_instance_f32 structure.
 */
void my_rfft_fast_init_2048(arm_rfft_fast_instance_f32 *S)
{
    // The CFFT instance is a nested struct 'Sint', not a pointer 'pCfft'.
    // We initialize it by copying the pre-defined constant struct.
    S->Sint = arm_cfft_sR_f32_len2048;
    S->fftLenRFFT = 2048;

    // The twiddle factors are now correctly declared via 'extern'.
    S->pTwiddleRFFT = (float32_t *) twiddleCoef_rfft_4096;
}

// ====================================================================
// 4. Core MFCC Processing for a single Frame
// ====================================================================

void ProcessFrame(float* current_frame) {
    // 1. Pre-emphasis
    float last_sample = prev_sample_for_preemphasis;
    for (int i = FRAME_LENGTH - 1; i > 0; --i) {
        current_frame[i] -= PRE_EMPHASIS_ALPHA * current_frame[i - 1];
    }
    current_frame[0] -= PRE_EMPHASIS_ALPHA * last_sample;
    prev_sample_for_preemphasis = current_frame[FRAME_LENGTH - 1]; // Save last sample for next frame

    // 2. Apply Hamming Window
    for (size_t j = 0; j < FRAME_LENGTH; ++j) {
        current_frame[j] *= hamming_window[j];
    }

    // 3. FFT (using CMSIS-DSP)
    // Zero-pad the frame into the FFT input buffer
    for(size_t j=0; j < FFT_SIZE; ++j) {
        fft_input[j] = (j < FRAME_LENGTH) ? current_frame[j] : 0.0f;
    }
    
    // Perform Real FFT
    arm_rfft_fast_f32(&fft_instance, fft_input, fft_output, 0);

    // 4. Power Spectrum
    // The output is interleaved [Re(0), Re(1), Im(1), ..., Re(N/2-1), Im(N/2-1), Re(N/2)]
    // We compute magnitude squared for each complex component
    arm_cmplx_mag_squared_f32(fft_output, power_spectrum, FFT_SIZE / 2 + 1);

    // 5. Apply Mel Filterbank
    for (int j = 0; j < NUM_MEL_BANDS; ++j) {
        mel_energies[j] = 0.0f;
        for (size_t k = 0; k < FFT_SIZE / 2 + 1; ++k) {
            mel_energies[j] += power_spectrum[k] * mel_filterbank[j][k];
        }
        // Apply log to the energies
        if (mel_energies[j] > 1e-6) {
             mel_energies[j] = logf(mel_energies[j]);
        } else {
             mel_energies[j] = -13.8155f; // log(1e-6)
        }
    }

    // 6. DCT-II
    dct_ii(mel_energies, mfcc_coeffs, NUM_MEL_BANDS, NUM_MFCC_COEFFS);
    
    // 7. Print the results
    hw.Print("MFCCs: [");
    for(int i = 0; i < NUM_MFCC_COEFFS; ++i) {
        hw.Print("%.3f, ", mfcc_coeffs[i]);
    }
    hw.PrintLine("]");
}


// ====================================================================
// 5. Audio Callback
// ====================================================================

void AudioCallback(AudioHandle::InputBuffer in, AudioHandle::OutputBuffer out, size_t size) {
    for (size_t i = 0; i < size; i++) {
        // We only care about the left input channel for this example
        circular_audio_buffer[write_pos] = in[0][i];
        
        // Advance write pointer
        write_pos = (write_pos + 1) % CIRCULAR_BUFFER_SIZE;
    }
}


// ====================================================================
// 6. Main Function
// ====================================================================

int main(void) {
    // Initialize the Daisy Seed hardware
    hw.Init();
    hw.SetAudioBlockSize(BLOCK_SIZE); 
    hw.SetAudioSampleRate(SaiHandle::Config::SampleRate::SAI_16KHZ);

    // Start serial logging
    hw.StartLog(true);
    hw.PrintLine("Starting Streaming MFCC Analyzer...");

    // Initialize MFCC components
    init_hamming_window();
    init_mel_filterbank();

    // Initialize the CMSIS-DSP FFT instance
    // if (arm_rfft_fast_init_f32(&fft_instance, FFT_SIZE) != ARM_MATH_SUCCESS) {
    //     hw.PrintLine("Error initializing FFT!");
    //     while(1) {}
    // }

    my_rfft_fast_init_2048(&fft_instance);

    
    // Start the audio callback
    hw.StartAudio(AudioCallback);

    while (1) {
        // Calculate how many new samples are available in the circular buffer
        int samples_available = (write_pos - read_pos + CIRCULAR_BUFFER_SIZE) % CIRCULAR_BUFFER_SIZE;

        // If we have enough samples to fill a frame, process it
        if (samples_available >= FRAME_LENGTH) {
            
            // Copy a frame's worth of data from the circular buffer
            for (int i = 0; i < FRAME_LENGTH; ++i) {
                frame[i] = circular_audio_buffer[(read_pos + i) % CIRCULAR_BUFFER_SIZE];
            }

            // Process this frame
            ProcessFrame(frame);

            // Advance the read pointer by the stride
            read_pos = (read_pos + FRAME_STRIDE) % CIRCULAR_BUFFER_SIZE;
        }

        // A small delay to prevent the main loop from running too fast
        System::Delay(1);
    }
}