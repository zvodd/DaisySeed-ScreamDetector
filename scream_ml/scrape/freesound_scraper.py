import scrapy
import argparse
from scrapy.crawler import CrawlerProcess
from scrapy.pipelines.files import FilesPipeline

# Define the item structure for the data to be scraped.
class MyFileItem(scrapy.Item):
    file_urls = scrapy.Field()
    files = scrapy.Field()

# The spider class defines how to crawl and extract data.
class FreesoundListSpider(scrapy.Spider):
    name = 'freesound_list_spider'

    def __init__(self, url_file=None, *args, **kwargs):
        """
        Initializes the spider, reading a list of URLs from the provided file path.
        """
        super(FreesoundListSpider, self).__init__(*args, **kwargs)
        self.start_urls = []
        if url_file:
            with open(url_file, 'r') as f:
                # Read each line, strip whitespace, and add non-empty lines to the list.
                self.start_urls = [line.strip() for line in f.readlines() if line.strip()]
        else:
            self.logger.error("No URL file provided! Use the 'url_file' argument.")

    def start_requests(self):
        """
        Generates the initial requests to be crawled from the start_urls list.
        """
        for url in self.start_urls:
            yield scrapy.Request(url=url, callback=self.parse)

    def parse(self, response):
        """
        Parses the response from each URL to find the audio file link.
        The link is located in the 'content' attribute of a <meta> tag.
        Example: <meta property="og:audio" content="https://cdn.freesound.org/.../file.mp3" />
        """
        # Use a direct CSS selector to find the content attribute of the specific meta tag.
        audio_url = response.css('meta[property="og:audio"]::attr(content)').get()

        if audio_url:
            if audio_url.startswith('https://freesound.orghttps://cdn.freesound.org'):
                audio_url = audio_url[len('https://freesound.org'):]
            self.logger.info(f'Found audio URL: {audio_url} on page {response.url}')
            
            # The FilesPipeline expects a list of URLs in the 'file_urls' field.
            item = MyFileItem()
            item['file_urls'] = [audio_url]
            yield item
        else:
            self.logger.warning(f'Could not find audio URL on page: {response.url}')


def run_spider(url_file_path):
    """
    Configures and runs the Scrapy spider process.
    """
    # Configure Scrapy settings. Most importantly, enable the FilesPipeline.
    settings = {
        'ITEM_PIPELINES': {'scrapy.pipelines.files.FilesPipeline': 1},
        'FILES_STORE': 'downloaded_audio',  # Specifies the directory to save files.
        'LOG_LEVEL': 'INFO', # Keeps the console output clean.
        'REQUEST_FINGERPRINTER_IMPLEMENTATION': '2.7', # Suppresses a deprecation warning
    }

    # Create a CrawlerProcess with the specified settings.
    process = CrawlerProcess(settings)

    # Add the spider to the process, passing the command-line argument to it.
    process.crawl(FreesoundListSpider, url_file=url_file_path)

    # Start the crawling process. The script will block here until done.
    print("--- Starting spider ---")
    process.start()
    print("--- Spider finished ---")


if __name__ == "__main__":
    # Set up command-line argument parsing.
    parser = argparse.ArgumentParser(
        description="A standalone Scrapy script to download audio files from a list of Freesound URLs."
    )
    # Add a required positional argument for the input file.
    parser.add_argument(
        "url_list_file",
        help="Path to the text file containing the URLs to download (one URL per line)."
    )
    args = parser.parse_args()

    # Run the spider with the provided file path.
    run_spider(url_file_path=args.url_list_file)