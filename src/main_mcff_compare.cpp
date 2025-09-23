#include "daisy_seed.h"
#include "arm_math.h"
#include <vector>

using namespace daisy;

DaisySeed hw;

// ====================================================================
// 0. About
// ====================================================================
// Pretty much a straight conversion of:
// https://www.geeksforgeeks.org/nlp/mel-frequency-cepstral-coefficients-mfcc-for-speech-recognition/

// There is an inital RMS threshold stage that triggers MFCC analisys to save power.

// So far the matching on MFCC-DCT-coefficients is untested and undecided.

// Each audio buffer would be 19.2 Kilobytes per 100ms of samples
// There are 2 300ms buffers, totaling ~120kb
// The rest of the buffers should be De minimis, say < 20kb





// ====================================================================
// 1. Configuration & Parameters (Unchanged)
// ====================================================================
const float SAMPLE_RATE = 48000.0f;
const size_t BLOCK_SIZE = 48;
const float RMS_THRESHOLD = 0.3f;
const uint32_t TRIGGER_HOLDOFF_MS = 1000;
const float BUFFER_DURATION_S = 0.3f;
const size_t BUFFER_SAMPLES = static_cast<size_t>(SAMPLE_RATE * BUFFER_DURATION_S);
const float PRE_EMPHASIS_ALPHA = 0.97f;
const size_t FRAME_LENGTH = static_cast<size_t>(0.025f * SAMPLE_RATE);
const size_t FRAME_STRIDE = static_cast<size_t>(0.010f * SAMPLE_RATE);
const size_t FFT_SIZE = 2048;
const int NUM_MEL_BANDS = 40;
const int NUM_MFCC_COEFFS = 13;
const size_t NUM_FRAMES = (BUFFER_SAMPLES > FRAME_LENGTH) ? (1 + (BUFFER_SAMPLES - FRAME_LENGTH) / FRAME_STRIDE) : 0;
// (BUFFER_SAMPLES[14400] - FRAME_LENGTH[1200]) / FRAME_STRIDE[480] = 27.5 Frames = 27


// ====================================================================
// 2. MFCC Pipeline Structure & Buffers
// ====================================================================

struct MfccPipeline
{
    // CMSIS-DSP instance for FFT
    arm_rfft_fast_instance_f32 fft_instance;

    // Pre-calculated tables
    float hamming_window[FRAME_LENGTH];
    std::vector<std::vector<float>> mel_filterbank;

    // Temporary processing buffers (kept in fast SRAM)
    float fft_input[FFT_SIZE];
    float fft_output[FFT_SIZE];
    float power_spectrum[FFT_SIZE / 2 + 1];
    float mel_energies[NUM_MEL_BANDS];
    float linear_audio_buffer[BUFFER_SAMPLES];
};

// ** Global instance of our pipeline
MfccPipeline mfcc;

// ** Buffers in large SDRAM
float DSY_SDRAM_BSS audio_buffer[BUFFER_SAMPLES];
float DSY_SDRAM_BSS calculated_mfccs[NUM_FRAMES][NUM_MFCC_COEFFS];

// Global state variables
size_t buffer_write_pos = 0;
volatile bool analysis_pending = false;
uint32_t last_trigger_time = 0;

// Reference MFCCs (placeholder)
const float reference_scream_mfccs[NUM_FRAMES][NUM_MFCC_COEFFS] = {{0.0f}};
const float SCREAM_MATCH_THRESHOLD = 10.0f;


// ====================================================================
// 3. Encapsulated MFCC Pipeline Functions
// ====================================================================

// Simple DCT-II function (as required)
void dct_ii(float* input, float* output, int size, int num_coeffs)
{
    for(int k = 0; k < num_coeffs; ++k)
    {
        float sum = 0.0f;
        for(int n = 0; n < size; ++n)
        {
            sum += input[n] * cosf(M_PI / size * (n + 0.5f) * k);
        }
        output[k] = sum;
    }
}

/**
 * @brief Initializes the MFCC pipeline (FFT, window, filterbank).
 * @param S Pointer to the MfccPipeline instance.
 */
void mfcc_pipeline_init(MfccPipeline* S)
{
    // Init CMSIS-DSP FFT
    arm_rfft_fast_init_f32(&S->fft_instance, FFT_SIZE);

    // Init Hamming Window
    for(size_t i = 0; i < FRAME_LENGTH; ++i)
    {
        S->hamming_window[i] = 0.54f - 0.46f * cosf(2.0f * M_PI * i / (FRAME_LENGTH - 1));
    }

    // Init Mel Filterbank (using helper functions)
    auto hz_to_mel = [](float hz) { return 2595.0f * log10f(1.0f + hz / 700.0f); };
    auto mel_to_hz = [](float mel) { return 700.0f * (powf(10.0f, mel / 2595.0f) - 1.0f); };

    S->mel_filterbank.resize(NUM_MEL_BANDS, std::vector<float>(FFT_SIZE / 2 + 1, 0.0f));
    float low_mel = hz_to_mel(0);
    float high_mel = hz_to_mel(SAMPLE_RATE / 2.0f);
    std::vector<float> mel_points(NUM_MEL_BANDS + 2);
    for(int i = 0; i < NUM_MEL_BANDS + 2; ++i)
    {
        mel_points[i] = low_mel + i * (high_mel - low_mel) / (NUM_MEL_BANDS + 1);
    }
    std::vector<size_t> fft_bins(NUM_MEL_BANDS + 2);
    for(int i = 0; i < NUM_MEL_BANDS + 2; ++i)
    {
        fft_bins[i] = static_cast<size_t>(floorf((FFT_SIZE + 1) * mel_to_hz(mel_points[i]) / SAMPLE_RATE));
    }
    for(int m = 1; m < NUM_MEL_BANDS + 1; ++m)
    {
        for(size_t k = fft_bins[m - 1]; k < fft_bins[m]; ++k)
        {
            S->mel_filterbank[m-1][k] = (float)(k - fft_bins[m-1]) / (fft_bins[m] - fft_bins[m-1]);
        }
        for(size_t k = fft_bins[m]; k < fft_bins[m+1]; ++k)
        {
            S->mel_filterbank[m-1][k] = (float)(fft_bins[m+1] - k) / (fft_bins[m+1] - fft_bins[m]);
        }
    }
}


/**
 * @brief Processes an audio buffer to produce MFCCs. This is our replacement.
 * @param S         Pointer to the initialized MfccPipeline instance.
 * @param audio_in  Pointer to the input audio buffer (size: BUFFER_SAMPLES).
 * @param mfcc_out  Pointer to the output 2D array for MFCC coefficients.
 */
void mfcc_pipeline_process(MfccPipeline* S, const float* audio_in, float (*mfcc_out)[NUM_MFCC_COEFFS])
{
    // 1. Copy audio and apply pre-emphasis
    arm_copy_f32(const_cast<float*>(audio_in), S->linear_audio_buffer, BUFFER_SAMPLES);
    for(int i = BUFFER_SAMPLES - 1; i > 0; --i)
    {
        S->linear_audio_buffer[i] -= PRE_EMPHASIS_ALPHA * S->linear_audio_buffer[i - 1];
    }

    // 2. Process each frame
    for(size_t i = 0; i < NUM_FRAMES; ++i)
    {
        size_t frame_start = i * FRAME_STRIDE;
        float* frame = &S->linear_audio_buffer[frame_start];

        // Apply Hamming Window
        arm_mult_f32(frame, S->hamming_window, frame, FRAME_LENGTH);

        // Prepare for FFT (copy and zero-pad)
        arm_fill_f32(0.0f, S->fft_input, FFT_SIZE);
        arm_copy_f32(frame, S->fft_input, FRAME_LENGTH);
        
        // 3. FFT
        arm_rfft_fast_f32(&S->fft_instance, S->fft_input, S->fft_output, 0);

        // 4. Power Spectrum
        arm_cmplx_mag_squared_f32(S->fft_output, S->power_spectrum, FFT_SIZE / 2 + 1);
        arm_scale_f32(S->power_spectrum, 1.0f / FFT_SIZE, S->power_spectrum, FFT_SIZE / 2 + 1);

        // 5. Apply Mel Filterbank & Log
        for(int j = 0; j < NUM_MEL_BANDS; ++j)
        {
            arm_dot_prod_f32(S->power_spectrum, S->mel_filterbank[j].data(), FFT_SIZE / 2 + 1, &S->mel_energies[j]);
            S->mel_energies[j] = (S->mel_energies[j] > 1e-6) ? logf(S->mel_energies[j]) : -13.8155f;
        }

        // 6. DCT
        dct_ii(S->mel_energies, mfcc_out[i], NUM_MEL_BANDS, NUM_MFCC_COEFFS);
    }
}


// ====================================================================
// 4. Main Logic
// ====================================================================

// Audio callback remains lean and fast
void audio_callback(AudioHandle::InputBuffer in, AudioHandle::OutputBuffer out, size_t size)
{
    uint32_t now = System::GetNow();
    arm_copy_f32(const_cast<float*>(in[0]), out[0], size);
    arm_copy_f32(const_cast<float*>(in[1]), out[1], size);

    for(size_t i = 0; i < size; ++i)
    {
        audio_buffer[buffer_write_pos] = in[0][i];
        buffer_write_pos = (buffer_write_pos + 1) % BUFFER_SAMPLES;
    }

    if(!analysis_pending && (now - last_trigger_time > TRIGGER_HOLDOFF_MS))
    {
        float rms_level = 0;
        arm_rms_f32(const_cast<float*>(in[0]), size, &rms_level);
        if(rms_level > RMS_THRESHOLD)
        {
            analysis_pending = true;
            last_trigger_time = now;
        }
    }
}

bool compare_mfccs()
{
    float total_error = 0.0f;
    for(size_t i = 0; i < NUM_FRAMES; ++i)
    {
        for(size_t j = 0; j < NUM_MFCC_COEFFS; ++j)
        {
            float diff = reference_scream_mfccs[i][j] - calculated_mfccs[i][j];
            total_error += diff * diff;
        }
    }
    float mse = total_error / (NUM_FRAMES * NUM_MFCC_COEFFS);

    return mse < SCREAM_MATCH_THRESHOLD;
}

// TODO: Exeriment with suitable matching techniques on series of Mel-DCTs
// void compare_delta_mfccs(){
//     // calculate the delta of consecutive frames, match when from
//     // delta0 = calculated_mfccs[1][:] - calculated_mfccs[0][:]

//     for(size_t i = 0; i < NUM_FRAMES; ++i)
//     {
//         // delta_w
//         for(size_t j = 0; j < NUM_MFCC_COEFFS; ++j)
//         {
//             //
//         }
//     }

// }

int main(void)
{
    hw.Init();
    hw.SetAudioBlockSize(BLOCK_SIZE);
    hw.SetAudioSampleRate(SaiHandle::Config::SampleRate::SAI_48KHZ);

    // ** Initialize the entire MFCC pipeline **
    mfcc_pipeline_init(&mfcc);

    hw.StartAudio(audio_callback);

    while(1)
    {
        if(analysis_pending)
        {
            // Create a linear copy of the circular buffer
            float temp_audio_buffer[BUFFER_SAMPLES]; 
            size_t read_pos = buffer_write_pos;
            for(int i = BUFFER_SAMPLES - 1; i >= 0; --i)
            {
                read_pos = (read_pos == 0) ? BUFFER_SAMPLES - 1 : read_pos - 1;
                temp_audio_buffer[i] = audio_buffer[read_pos];
            }
            
            // ** Execute the entire pipeline with a single function call **
            mfcc_pipeline_process(&mfcc, temp_audio_buffer, calculated_mfccs);
            
            // Compare the result
            if (compare_mfccs())
            {
                // detection event
                hw.SetLed(true);
            }
            else
            {
                hw.SetLed(false);
            }
            // Reset the flag
            analysis_pending = false;
        }
        
        System::Delay(5);
    }
}