import numpy as np
import librosa
import librosa.display
import matplotlib.pyplot as plt
from matplotlib.widgets import Cursor
from scipy.fftpack import dct
from abc import ABC, abstractmethod


class AudioAnalyzer(ABC):
    """Abstract base class for audio analysis methods"""
    
    @abstractmethod
    def analyze(self, y, sr):
        """
        Analyze audio signal and return 2D array
        Args:
            y: audio time series
            sr: sample rate
        Returns:
            2D numpy array (features x time_frames)
        """
        pass
    
    @abstractmethod
    def get_feature_labels(self):
        """Return list of feature labels for y-axis"""
        pass
    
    @abstractmethod
    def get_title(self):
        """Return title for the plot"""
        pass


class MFCCAnalyzer(AudioAnalyzer):
    """MFCC analysis implementation"""
    
    def __init__(self, n_mfcc=13, hop_length=512, n_fft=2048, pre_emphasis=0.97):
        self.n_mfcc = n_mfcc
        self.hop_length = hop_length
        self.n_fft = n_fft
        self.pre_emphasis = pre_emphasis
    
    def analyze(self, y, sr):
        # Your existing MFCC computation or use librosa
        mfcc = librosa.feature.mfcc(y=y, sr=sr, n_mfcc=self.n_mfcc, 
                                   hop_length=self.hop_length, n_fft=self.n_fft)
        return mfcc
    
    def get_feature_labels(self):
        return [f'MFCC {i}' for i in range(self.n_mfcc)]
    
    def get_title(self):
        return 'MFCC Features'


class SpectrogramAnalyzer(AudioAnalyzer):
    """Spectrogram analysis implementation"""
    
    def __init__(self, hop_length=512, n_fft=2048):
        self.hop_length = hop_length
        self.n_fft = n_fft
    
    def analyze(self, y, sr):
        # Compute magnitude spectrogram
        stft = librosa.stft(y, hop_length=self.hop_length, n_fft=self.n_fft)
        magnitude = np.abs(stft)
        # Convert to dB
        magnitude_db = librosa.amplitude_to_db(magnitude, ref=np.max)
        return magnitude_db
    
    def get_feature_labels(self):
        return [f'Freq {i}' for i in range(self.n_fft // 2 + 1)]
    
    def get_title(self):
        return 'Spectrogram (dB)'


class MelSpectrogramAnalyzer(AudioAnalyzer):
    """Mel-spectrogram analysis implementation"""
    
    def __init__(self, n_mels=128, hop_length=512, n_fft=2048):
        self.n_mels = n_mels
        self.hop_length = hop_length
        self.n_fft = n_fft
    
    def analyze(self, y, sr):
        mel_spec = librosa.feature.melspectrogram(y=y, sr=sr, n_mels=self.n_mels,
                                                 hop_length=self.hop_length, n_fft=self.n_fft)
        mel_spec_db = librosa.power_to_db(mel_spec, ref=np.max)
        return mel_spec_db
    
    def get_feature_labels(self):
        return [f'Mel {i}' for i in range(self.n_mels)]
    
    def get_title(self):
        return 'Mel Spectrogram (dB)'


class AudioAnalysisPipeline:
    """Main pipeline class for audio analysis and visualization"""
    
    def __init__(self, analyzer: AudioAnalyzer):
        self.analyzer = analyzer
        self.audio_data = None
        self.sr = None
        self.analysis_result = None
        self.selected_frame = None
        self.selected_values = None
        
        # Visualization attributes
        self.fig = None
        self.ax_wave = None
        self.ax_analysis = None
        self.cursor = None
        self.vertical_line = None
        
    def load_audio(self, file_path, sr=None):
        """Load audio file"""
        self.audio_data, self.sr = librosa.load(file_path, sr=sr)
        print(f"Loaded audio: {len(self.audio_data)} samples at {self.sr} Hz")
        return self.audio_data, self.sr
    
    def analyze(self):
        """Run analysis on loaded audio"""
        if self.audio_data is None:
            raise ValueError("No audio data loaded. Call load_audio() first.")
        
        self.analysis_result = self.analyzer.analyze(self.audio_data, self.sr)
        print(f"Analysis result shape: {self.analysis_result.shape}")
        return self.analysis_result
    
    def visualize(self, figsize=(14, 8)):
        """Create interactive visualization"""
        if self.analysis_result is None:
            raise ValueError("No analysis results. Call analyze() first.")
        
        self.fig, (self.ax_wave, self.ax_analysis) = plt.subplots(2, 1, figsize=figsize)
        
        # Plot waveform
        librosa.display.waveshow(self.audio_data, sr=self.sr, ax=self.ax_wave)
        self.ax_wave.set_title('Waveform')
        self.ax_wave.set_xlabel('Time (s)')
        
        # Plot analysis result
        im = self.ax_analysis.imshow(self.analysis_result, cmap='hot', aspect='auto', 
                                   origin='lower', interpolation='nearest')
        self.ax_analysis.set_title(self.analyzer.get_title())
        self.ax_analysis.set_xlabel('Frame Index')
        self.ax_analysis.set_ylabel('Feature Index')
        
        # Add colorbar
        plt.colorbar(im, ax=self.ax_analysis)
        
        # Add cursor for interaction
        self.cursor = Cursor(self.ax_analysis, useblit=True, color='white', linewidth=2)
        
        # Connect click event
        self.fig.canvas.mpl_connect('button_press_event', self._on_click)
        
        plt.tight_layout()
        plt.show()
        
        return self.fig
    
    def _on_click(self, event):
        """Handle click events on the analysis plot"""
        if event.inaxes == self.ax_analysis:
            frame_idx = int(round(event.xdata))
            
            if 0 <= frame_idx < self.analysis_result.shape[1]:
                self.selected_frame = frame_idx
                self.selected_values = self.analysis_result[:, frame_idx].copy()
                
                # Remove previous vertical line
                if self.vertical_line is not None:
                    self.vertical_line.remove()
                
                # Add new vertical line
                self.vertical_line = self.ax_analysis.axvline(x=frame_idx, color='cyan', 
                                                            linewidth=2, alpha=0.8)
                
                # Update display
                self.fig.canvas.draw()
                
                print(f"\nSelected frame: {frame_idx}")
                print(f"Feature values shape: {self.selected_values.shape}")
                print("Feature values:")
                feature_labels = self.analyzer.get_feature_labels()
                for i, value in enumerate(self.selected_values):
                    label = feature_labels[i] if i < len(feature_labels) else f"Feature {i}"
                    print(f"  {label}: {value:.4f}")
    
    def get_selected_values(self):
        """Get the currently selected vertical slice values"""
        if self.selected_values is None:
            print("No frame selected. Click on the analysis plot to select a frame.")
            return None
        return self.selected_frame, self.selected_values.copy()
    
    def plot_selected_values(self):
        """Plot the selected vertical slice as a separate figure"""
        if self.selected_values is None:
            print("No frame selected. Click on the analysis plot to select a frame.")
            return None
        
        fig, ax = plt.subplots(figsize=(10, 6))
        feature_labels = self.analyzer.get_feature_labels()
        
        x_positions = range(len(self.selected_values))
        ax.plot(x_positions, self.selected_values, 'o-', linewidth=2, markersize=6)
        
        ax.set_xlabel('Feature Index')
        ax.set_ylabel('Feature Value')
        ax.set_title(f'Selected Features - Frame {self.selected_frame}')
        ax.grid(True, alpha=0.3)
        
        # Add feature labels if not too many
        if len(self.selected_values) <= 20:
            ax.set_xticks(x_positions[::max(1, len(x_positions)//10)])
            ax.set_xticklabels([feature_labels[i] if i < len(feature_labels) 
                               else f"F{i}" for i in x_positions[::max(1, len(x_positions)//10)]], 
                              rotation=45)
        
        plt.tight_layout()
        plt.show()
        return fig


# Example usage and demonstration
if __name__ == "__main__":
    # Example of how to use the pipeline
    
    # Create different analyzers
    mfcc_analyzer = MFCCAnalyzer(n_mfcc=13)
    # mel_analyzer = MelSpectrogramAnalyzer(n_mels=40)
    # spec_analyzer = SpectrogramAnalyzer()
    
    # # Example usage (commented out since we don't have an actual audio file)
    # """
    # Create pipeline with MFCC analyzer
    pipeline = AudioAnalysisPipeline(mfcc_analyzer)
    
    # Load audio file
    audio_data, sr = pipeline.load_audio('./freesounds/608752_13454867-hq.mp3')
    
    # Analyze
    analysis_result = pipeline.analyze()
    
    # Visualize (creates interactive plot)
    pipeline.visualize()
    
    # After clicking on the plot, get selected values
    frame_idx, selected_values = pipeline.get_selected_values()
    
    # Plot the selected slice
    pipeline.plot_selected_values()
    
    # Switch to different analyzer
    pipeline.analyzer = mel_analyzer
    pipeline.analyze()  # Re-analyze with new method
    pipeline.visualize()  # New interactive plot
    # """
    
    # print("Pipeline classes created successfully!")
    # print("\nTo use:")
    # print("1. Create an analyzer: analyzer = MFCCAnalyzer()")
    # print("2. Create pipeline: pipeline = AudioAnalysisPipeline(analyzer)")
    # print("3. Load audio: pipeline.load_audio('file.wav')")
    # print("4. Analyze: pipeline.analyze()")
    # print("5. Visualize: pipeline.visualize()")
    # print("6. Click on the plot to select a frame")
    # print("7. Get values: frame, values = pipeline.get_selected_values()")