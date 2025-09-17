#include "daisy_seed.h"
#include <math.h>

using namespace daisy;

DaisySeed hw;

// Audio processing parameters
const float SAMPLE_RATE = 48000.0f;
const float THRESHOLD = 0.3f;  // RMS threshold (0.0 to 1.0)
const size_t BLOCK_SIZE = 48;  // Audio block size

// LED pulse parameters  
const uint32_t PULSE_DURATION_MS = 100;
uint32_t pulse_start_time = 0;
bool led_active = false;

// RMS calculation variables
float rms_buffer[BLOCK_SIZE];
size_t buffer_index = 0;

// Simple RMS calculation over audio block
float calculate_rms(float* audio_block, size_t block_size) {
    float sum_squares = 0.0f;
    
    for(size_t i = 0; i < block_size; i++) {
        float sample = audio_block[i];
        sum_squares += sample * sample;
    }
    
    return sqrtf(sum_squares / block_size);
}

// Audio callback function
void audio_callback(AudioHandle::InputBuffer in, 
                   AudioHandle::OutputBuffer out, 
                   size_t size) {
    
    // Process left channel input
    const float* input = in[0];
    
    // Calculate RMS of current audio block
    float rms_level = calculate_rms(const_cast<float*>(input), size);
    
    // Check if RMS exceeds threshold
    if(rms_level > THRESHOLD && !led_active) {
        // Trigger LED pulse
        led_active = true;
        pulse_start_time = System::GetNow();
        hw.SetLed(1); // Turn LED on
    }
    
    // Pass audio through (optional - remove if no audio output needed)
    for(size_t i = 0; i < size; i++) {
        out[0][i] = in[0][i]; // Left channel pass-through
        out[1][i] = in[1][i]; // Right channel pass-through
    }
}

int main(void) {
    hw.Configure();
    // Initialize hardware
    hw.Init();
    hw.SetAudioBlockSize(BLOCK_SIZE);

    hw.SetAudioSampleRate(daisy::SaiHandle::Config::SampleRate::SAI_48KHZ);
    
    hw.SetLed(0);
    
    // Start audio processing
    hw.StartAudio(audio_callback);
    
    // Main loop
    while(1) {
        // Handle LED pulse timing
        if(led_active) {
            uint32_t current_time = System::GetNow();
            if(current_time - pulse_start_time >= PULSE_DURATION_MS) {
                hw.SetLed(0); // Turn LED off
                led_active = false;
            }
        }
        
        
        // Small delay to prevent excessive CPU usage
        System::Delay(1);
    }
}