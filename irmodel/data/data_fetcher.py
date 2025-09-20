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
        """
        Fetches the zero-coupon bond price data from the RDA file.
        Returns:
            pd.DataFrame: DataFrame containing:
                Time to Maturity (float64): Time to maturity in years.
                Strip (float64): Price of the zero-coupon bond.
        """
        result = pyreadr.read_r(self.data_dir + '/VeronesiTable15p1.rda')
        return result['VeronesiTable15p1']

    def hull_white_data(self) -> 'pd.DataFrame':
        """
        Fetches the cap price data from the RDA file.
        Returns:
            pd.DataFrame: DataFrame containing:
                Maturity (float64): Maturity in years.
                Swap Rate (float64): Strike price.
                Cap Price (x100) (float64): Price of the cap.
                Discount Factor (float64): Discount factor.
        """
        result = pyreadr.read_r(self.data_dir + '/VeronesiTable19p4.rda')
        return result['VeronesiTable19p4']


