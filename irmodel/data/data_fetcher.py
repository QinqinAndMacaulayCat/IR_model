import pyreadr
import os

class SampleDataFetcher:
    
    def __init__(self):
        self.data_dir = os.path.dirname(__file__)

    def short_rate(self) -> 'pd.DataFrame':
        """Fetches the short rate data from the RDA file.
        Returns:
            pd.DataFrame: DataFrame containing:
                - DATE (datetime64[ns]): Date of the observation.
                - rt (float64): Short rate value.
        """
        result = pyreadr.read_r(self.data_dir + '/VeronesiTable14p7q5.rda')  
        return result['VeronesiTable14p7q5']
    
    def zcb_price(self) -> 'pd.DataFrame':
        result = pyreadr.read_r(self.data_dir + '/VeronesiTable15p1.rda')
        return result['VeronesiTable15p1']


if __name__ == "__main__":
    fetcher = SampleDataFetcher()
    print(fetcher.short_rate())
    print(fetcher.short_rate().dtypes)
    print(fetcher.zcb_price())

