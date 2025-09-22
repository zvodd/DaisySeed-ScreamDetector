#include "daisy_seed.h"
#include <math.h>

using namespace daisy;

DaisySeed hw;
Led led;

// Audio processing parameters
const float SAMPLE_RATE = 48000.0f;
const float THRESHOLD = 0.1f;  // RMS threshold (0.0 to 1.0)
const size_t BLOCK_SIZE = 48;  // Audio block size (1ms at 48kHz)

// Analysis window parameters
const size_t ANALYSIS_WINDOW_MS = 400;  // 400ms analysis window
const size_t BLOCKS_PER_WINDOW = (ANALYSIS_WINDOW_MS * SAMPLE_RATE / 1000) / BLOCK_SIZE; // 400 blocks
const size_t TOTAL_SAMPLES_IN_WINDOW = BLOCKS_PER_WINDOW * BLOCK_SIZE; // 19,200 samples

// LED pulse parameters  
const uint32_t PULSE_DURATION_MS = 100;
uint32_t pulse_start_time = 0;
bool led_active = false;

// Circular buffer for streaming audio analysis
class AudioCircularBuffer {
private:
    float* buffer;
    size_t buffer_size;
    size_t write_index;
    size_t blocks_written;
    
public:
    AudioCircularBuffer(size_t size) : buffer_size(size), write_index(0), blocks_written(0) {
        buffer = new float[buffer_size];
        // Initialize with zeros
        for(size_t i = 0; i < buffer_size; i++) {
            buffer[i] = 0.0f;
        }
    }
    
    ~AudioCircularBuffer() {
        delete[] buffer;
    }
    
    // Add a new audio block to the circular buffer
    void add_block(const float* audio_block, size_t block_size) {
        for(size_t i = 0; i < block_size; i++) {
            buffer[write_index] = audio_block[i];
            write_index = (write_index + 1) % buffer_size;
        }
        blocks_written++;
    }
    
    // Check if we have enough data for analysis (full window)
    bool is_window_ready() const {
        return blocks_written >= BLOCKS_PER_WINDOW;
    }
    
    // Get the current analysis window (400ms worth of samples)
    // Returns samples in chronological order (oldest to newest)
    void get_analysis_window(float* output_window) const {
        if (!is_window_ready()) return;
        
        size_t read_start = write_index; // Start from oldest sample
        
        for(size_t i = 0; i < buffer_size; i++) {
            size_t read_index = (read_start + i) % buffer_size;
            output_window[i] = buffer[read_index];
        }
    }
    
    // Get recent RMS over last N blocks for quick threshold checking
    float get_recent_rms(size_t num_blocks = 10) const {
        if (blocks_written == 0) return 0.0f;
        
        size_t samples_to_check = num_blocks * BLOCK_SIZE;
        if (samples_to_check > buffer_size) samples_to_check = buffer_size;
        
        float sum_squares = 0.0f;
        size_t start_index = (write_index - samples_to_check + buffer_size) % buffer_size;
        
        for(size_t i = 0; i < samples_to_check; i++) {
            size_t idx = (start_index + i) % buffer_size;
            float sample = buffer[idx];
            sum_squares += sample * sample;
        }
        
        return sqrtf(sum_squares / samples_to_check);
    }
    
    size_t get_blocks_written() const { return blocks_written; }
    size_t get_buffer_size() const { return buffer_size; }
};

// Global circular buffer instance
AudioCircularBuffer* audio_buffer = nullptr;

// Analysis window for MFCC or other feature extraction
float* analysis_window = nullptr;

// Simple RMS calculation over audio block (for immediate threshold checking)
float calculate_rms(const float* audio_block, size_t block_size) {
    float sum_squares = 0.0f;
    
    for(size_t i = 0; i < block_size; i++) {
        float sample = audio_block[i];
        sum_squares += sample * sample;
    }
    
    return sqrtf(sum_squares / block_size);
}

// Placeholder for future MFCC or feature analysis
void analyze_audio_window(const float* window, size_t window_size) {
    // This is where you would implement:
    // - MFCC feature extraction
    // - Spectral analysis (FFT)
    // - Pattern recognition
    // - Machine learning inference
    // etc.
    
    // For now, just calculate overall RMS of the entire window
    float window_rms = calculate_rms(window, window_size);
    
    // Example: Could trigger different behaviors based on longer-term analysis
    if (window_rms > THRESHOLD * 0.5f) {
        // Detected sustained loud audio over 400ms window
        // Could implement different LED patterns, logging, etc.
    }
}

// Audio callback function
void audio_callback(AudioHandle::InputBuffer in, 
                   AudioHandle::OutputBuffer out, 
                   size_t size) {
    
    // Process left channel input
    const float* input = in[0];
    
    // Add current block to circular buffer
    audio_buffer->add_block(input, size);
    
    // Quick RMS check for immediate LED triggering (low latency)
    float immediate_rms = calculate_rms(input, size);
    
    // Check if immediate RMS exceeds threshold for LED trigger
    if(immediate_rms > THRESHOLD && !led_active) {
        led_active = true;
        pulse_start_time = System::GetNow();
        led.Set(1.0f);
    }
    
    // Perform larger window analysis if we have enough data
    if(audio_buffer->is_window_ready()) {
        // Every 10 blocks (10ms), perform window analysis
        static size_t analysis_counter = 0;
        analysis_counter++;
        
        if(analysis_counter >= 10) {
            analysis_counter = 0;
            
            // Get 400ms window for analysis
            audio_buffer->get_analysis_window(analysis_window);
            
            // Perform feature analysis (MFCC, etc.)
            analyze_audio_window(analysis_window, TOTAL_SAMPLES_IN_WINDOW);
        }
    }
    
    // Pass audio through (optional)
    for(size_t i = 0; i < size; i++) {
        out[0][i] = in[0][i]; // Left channel pass-through
        out[1][i] = in[1][i]; // Right channel pass-through
    }
}

int main(void) {
    // Initialize hardware
    hw.Init();
    hw.SetAudioBlockSize(BLOCK_SIZE);
    hw.SetAudioSampleRate(DaisySeed::AudioSampleRate::SAI_48KHZ);
    
    // Initialize LED
    led.Init(hw.GetPin(22), false);
    led.Set(0.0f);
    
    // Initialize circular buffer for 400ms of audio
    audio_buffer = new AudioCircularBuffer(TOTAL_SAMPLES_IN_WINDOW);
    
    // Allocate analysis window buffer
    analysis_window = new float[TOTAL_SAMPLES_IN_WINDOW];
    
    // Start audio processing
    hw.StartAudio(audio_callback);
    
    // Main loop
    while(1) {
        // Handle LED pulse timing
        if(led_active) {
            uint32_t current_time = System::GetNow();
            if(current_time - pulse_start_time >= PULSE_DURATION_MS) {
                led.Set(0.0f);
                led_active = false;
            }
        }
        
        // Update LED state
        led.Update();
        
        // Optional: Print buffer status for debugging
        static uint32_t last_debug_time = 0;
        uint32_t current_time = System::GetNow();
        if(current_time - last_debug_time > 1000) { // Every 1 second
            last_debug_time = current_time;
            // Could print buffer fill status, RMS levels, etc.
        }
        
        System::Delay(1);
    }
    
    // Cleanup (never reached in this application)
    delete audio_buffer;
    delete[] analysis_window;
}