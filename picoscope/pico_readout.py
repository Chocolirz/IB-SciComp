from picosdk.discover import find_unit
from picosdk.device import ChannelConfig, TimebaseOptions

__all__ = ['PicoReadout']

class PicoReadout:
    def __init__(self, ADecouple, BDecouple, ARange, BRange, min_collection_time) -> None:
        self.ADecouple = ADecouple
        self.BDecouple = BDecouple
        self.ARange = ARange
        self.BRange = BRange
        self.min_collection_time = min_collection_time

    def run(self, ADecouple, BDecouple, ARange, BRange, min_collection_time):
        with find_unit() as device:
            channel_configs = [ChannelConfig('A', True, ADecouple, ARange), 
                            ChannelConfig('B', True, BDecouple, BRange)]
            timebase_options = TimebaseOptions(min_collection_time = min_collection_time)
            times, voltages, overflow_warnings = device.capture_block(timebase_options, channel_configs)
        return times, voltages, overflow_warnings

    def __str__(self) -> str:
        return f'Picoscope Setting: ADecouple={self.ADecouple}, BDecouple={self.BDecouple}, ARange={self.ARange}, BRange={self.BRange}, min_collection_time={self.min_collection_time}'