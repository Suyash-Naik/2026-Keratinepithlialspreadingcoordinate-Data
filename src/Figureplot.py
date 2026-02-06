#import required libraries
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import os
from glob import glob
import pandas as pd
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D
from re import findall as find
import re
from datetime import datetime
import urllib.request
import urllib.error
today=datetime.today().strftime('%Y%m%d')
print("Figureplot.py run on ", today)
class Figureplot:
    def __init__(self, data_dir='data', output_dir='figures'):
        self.data_dir = data_dir
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)

    def download_file(self, url, figure_label):
        """
        Download a file from an HTTPS link and save it with a figure label.
        
        Parameters:
        -----------
        url : str
            The HTTPS URL of the file to download
        figure_label : str
            A label for the file (used in the filename)
            
        Returns:
        --------
        str
            The path to the downloaded file
            
        Raises:
        -------
        urllib.error.URLError
            If the download fails
        """
        try:
            # Extract file extension from URL
            file_extension = url.split('.')[-1].split('?')[0]
            
            # Create filename with figure label and timestamp
            filename = f"{figure_label}_{today}.{file_extension}"
            filepath = os.path.join(self.data_dir, filename)
            
            # Create data directory if it doesn't exist
            os.makedirs(self.data_dir, exist_ok=True)
            
            # Download the file
            print(f"Downloading from: {url}")
            urllib.request.urlretrieve(url, filepath)
            
            print(f"File downloaded successfully to: {filepath}")
            return filepath
            
        except urllib.error.URLError as e:
            print(f"Error downloading file from {url}: {e}")
            raise
        except Exception as e:
            print(f"Unexpected error: {e}")
            raise