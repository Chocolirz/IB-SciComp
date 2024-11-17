###########
# This code was written during July & August 2022 by sfb47, drawing heavily for guidance on a piece of code
# ('EIS_asa73') by Arun Atwal, https://github.com/Arun-Atwal. It works on the Part IB picoscope (model 2206B) however
# it does not work on the Part II (2205A) or IA (3204) models because of an issue with picosdk v1.1 which prevents
# them being opened, and the IA scope in any case has poorly defined voltages ranges which raise a zero-indexing
# error. There also exists a program (Manual_PSD.ipynb) using the same demodulating function which requires manual
# upload of data.
###########

# Import Modules
import sys
from time import sleep
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.animation as animation
import tkinter as tk
from tkinter import messagebox
import picosdk
from picosdk import discover
from picosdk import device
import math

# INITIALISE THE GUI
class gui(tk.Tk):
    def __init__(self):
        tk.Tk.__init__(self)
        self._frame = None

        global defaults, scope
        scope = None
        defaults = True
        try:
            scope = discover.find_unit()
        except:
            pass
        if scope == None:
            self.switch_frame(NoDevicePage)
        else:
            self.switch_frame(StartPage)

    def switch_frame(self, frame_class):    # Function to switch between frames
        new_frame = frame_class(self)
        if self._frame is not None:
            self._frame.destroy()
        self._frame = new_frame
        self._frame.pack()

# PAGE WHEN NO DEVICE FOUND
class NoDevicePage(tk.Frame):
    def __init__(self, master):
        tk.Frame.__init__(self, master)
        self.errorlbl = tk.Label(self, text="No device found", font="Arial 20", fg="red").pack()
        tk.Button(self, text="Search for device", command=lambda: self.retry_press(master)).pack()

    def retry_press(self, master):  # Button to try again to find device
        global scope
        try:
            scope = discover.find_unit()
        except:
            pass
        if scope is not None:
            master.switch_frame(StartPage)

class Scope_Parameters():   # Pretty self-explanatory what this class is for.
                            # Note that it does not have the same order as the parameter entry page
    def __init__(self, mod_channel, VrangeA, VrangeB, f, totcyc, avgcyc, A, B):
        self.modulation_channel = mod_channel   # Channel containing modulating signal (string: A, B)
        self.Vinput_rangeA = VrangeA            # +/-V input range for channel A (float: from set of Vinput_ranges)
        self.Vinput_rangeB = VrangeB            # +/-V input range for channel B (float: from set of Vinput_ranges)
        self.frequency = f                      # Modulating frequency (float)
        self.total_cycles = totcyc              # Total number of cycles to record (float)
        self.averaging_cycles = avgcyc          # Number of cycles to average over (integer)
        self.Amode = A                          # Channel A mode (string: AC, DC)
        self.Bmode = B                          # Channel A mode (string: AC, DC)

        self.period = 1/self.frequency
        self.averaging_time = self.period * self.averaging_cycles
        self.sampling_time = self.period * self.total_cycles
        self.Nsamples = 5000


# PAGE TO SET SCOPE PARAMETERS
class StartPage(tk.Frame):
    def __init__(self, master):
        tk.Frame.__init__(self, master)

        # Read out device ID
        self.deviceID = tk.Frame()
        self.deviceID.pack()

        global scope
        #device_series_raw = str(scope.info.driver)
        device_model_raw = str(scope.info.variant)
        #device_series = device_series_raw.replace('picosdk ', '').replace('library', 'series;')
        device_model = device_model_raw.replace("'", '').replace('b', '')
        tk.Label(master=self.deviceID, text=f"Device found: PicoScope {device_model}", font="Arial 10", fg="black").pack()

        # CREATE PARAMETER-SETTING BOXES
        self.frm = tk.Frame()
        self.frm.pack()
        self.entries = []     # Empty array of parameter boxes
        startpagefont = 'Arial 14'

        # Labels to appear next to input boxes
        self.scope_parameters_names = [ 'Modulation signal channel:',
                                        'Channel A voltage range (\u00B1):',
                                        'A channel mode:',
                                        'Channel B voltage range (\u00B1):',
                                        'B channel mode:',
                                        'Modulating frequency (Hz):',
                                        'Total cycles:',
                                        'Averaging cycles:']

        # Names of input variables
        self.scope_parameters = ['mod_channel',
                                 'VrangeA',
                                 'Amode',
                                 'VrangeB',
                                 'Bmode',
                                 'f',
                                 'totcycles',
                                 'avgcycles']

        # Options available for inputs
        self.scope_parameters_options = [['A', 'B'],
                                         ['20mV', '50mV', '100mV', '200mV', '500mV', '1V', '2V', '5V', '10V', '20V'],
                                         ['DC', 'AC'],
                                         ['20mV', '50mV', '100mV', '200mV', '500mV', '1V', '2V', '5V', '10V', '20V'],
                                         ['DC', 'AC'],
                                         [],
                                         [],
                                         []]

        global defaults
        # Default entries
        self.default_parameters  = ['A',   # modulating channel
                                    '5V', # A voltage range
                                    'DC',  # channel A mode
                                    '100mV', # B voltage range
                                    'AC',  # channel B mode
                                    '',    # frequency
                                    '60', # total cycles
                                    '20']  # averaging cycles
        if defaults:    # Reset to the above values
            self.default_scope_parameters = self.default_parameters
        else:           # Keep the same values as already input
            Vinput_ranges_keys = {0.02: '20mV', 0.05: '50mV', 0.1: '100mV', 0.2: '200mV', 0.5: '500mV',
                                  1.0: '1V', 2.0: '2V', 5.0: '5V', 10.0: '10V', 20.0: '20V'}
            self.default_scope_parameters = [params.modulation_channel,
                                             Vinput_ranges_keys[params.Vinput_rangeA],
                                             params.Amode,
                                             Vinput_ranges_keys[params.Vinput_rangeB],
                                             params.Bmode,
                                             params.frequency,
                                             params.total_cycles,
                                             params.averaging_cycles]

        # Create the parameter-entry labels and buttons
        for i in range (0,8):
            label = tk.Label(master=self.frm, text=self.scope_parameters_names[i], font=startpagefont, fg='black',
                             anchor=tk.W, justify=tk.LEFT)  # Create input label
            if i in [5,6,7]:                                # Entry-type variables
                entry = tk.Entry(master=self.frm, font='Arial 12', justify=tk.RIGHT, width=10)  # Create entry box
                self.entries.append(entry)                          # Add entry to the array
                entry.insert(0, self.default_scope_parameters[i])   # Set default value
            else:   # Menu-type variables
                self.scope_parameters[i] = tk.StringVar(self)       # Initialise variable
                self.scope_parameters[i].set(self.default_scope_parameters[i])  # Set default value
                entry = tk.OptionMenu(self.frm, self.scope_parameters[i], *self.scope_parameters_options[i])  # Create dropdown box
                self.entries.append(entry)
            label.grid(row=len(self.entries), column=0, sticky=tk.W, pady=3)  # Place the label
            entry.grid(row=len(self.entries), column=1, sticky=tk.E, pady=3)  # Place the entry

        # Create buttons at bottom of page
        self.frm_buttons = tk.Frame()
        self.frm_buttons.pack(fill=tk.X, padx=5, pady=5)

        # Reset parameters to default
        btn_reset = tk.Button(master=self.frm_buttons, text = 'Reset to defaults', command = self.reset_parameters)
        btn_reset.pack(side=tk.LEFT, padx=10)
        # Start scope view
        btn_start = tk.Button(master=self.frm_buttons, text = 'Start', font = 'bold', bg = 'green', command = lambda: self.begin(master))
        btn_start.pack(side=tk.RIGHT, padx=10)

    # FUNCTIONS DEFINING THE ACTIONS OF THE BUTTONS
    def reset_parameters(self): # Resets parameters to defaults
        self.default_scope_parameters = self.default_parameters
        for i in range(0,8):
            if i in [5,6,7]:         # Return the entry-type options to their defaults
                self.entries[i].delete(0,'end')
                self.entries[i].insert(0, self.default_scope_parameters[i])
            else:   # Return the menu-type options to their defaults
                self.scope_parameters[i].set(self.default_scope_parameters[i])

    def check_params(self, params): # Check that the inputs are reasonable
        if params.frequency <= 0:
            return False, 'Frequency must be positive'
        if params.averaging_cycles <= 0:
            return False, 'Number of cycles must be positive'
        if params.total_cycles <= 0:
            return False, 'Number of cycles must be positive'
        if params.averaging_cycles > params.total_cycles:
            return False, 'Averaging cycles must be less than total'
        # If there are no problems:
        return True, None

    def begin(self, master):    # Switch to the scope view (get cracking)
        # Save the desired scope parameters as an object of type Scope_Parameters
        Vinput_ranges_floats = {'20mV': 20e-3, '50mV': 50e-3, '100mV': 100e-3, '200mV': 200e-3, '500mV': 500e-3,
                                '1V': 1., '2V': 2., '5V': 5., '10V': 10., '20V': 20.}  # Dictionary of allowed voltage ranges as floats
        entry_types_valid = False
        try:    # Set the scope parameters from the inputs
            global params
            params = Scope_Parameters(
                mod_channel=    str(self.scope_parameters[0].get()),
                VrangeA=        float(Vinput_ranges_floats[self.scope_parameters[1].get()]),
                VrangeB=        float(Vinput_ranges_floats[self.scope_parameters[3].get()]),
                f=              float(self.entries[5].get()),
                totcyc=         float(self.entries[6].get()),
                avgcyc=         int(self.entries[7].get()),
                A=              str(self.scope_parameters[2].get()),
                B=              str(self.scope_parameters[4].get())
            )
            entry_types_valid = True
        except:
            messagebox.showwarning(title="Input Error", message="Input(s) invalid.\nPlease check inputs!")

        if entry_types_valid:   # If inputs are reasonable
            entries_types_valid, fault = self.check_params(params)
            if entry_types_valid:   # If there are no problems
                # Destroy frames and open scope view
                self.deviceID.destroy()
                self.frm.destroy()
                self.frm_buttons.destroy()
                master.switch_frame(ScopePage)

            else:   # If inputs are not reasonable
                messagebox.showwarning(title="Input Error", message=f"Input(s) invalid.\nFault: {fault}\n")


# PAGE WITH THE SCOPE VIEWS
class ScopePage(tk.Frame):
    def __init__(self, master):
        tk.Frame.__init__(self, master)

        # Plots frame
        self.plots_frm = tk.Frame(self)
        self.plots_frm.grid(row = 0, column = 0)
        self.setup_plot(self.plots_frm)
        self.paused = False

        # Buttons frame
        self.button_frm = tk.Frame(self)
        self.button_frm.grid(row = 1, column = 0)
        return_btn = tk.Button(master=self.button_frm, text='Change parameters', font = 'Arial 10', fg = 'white', bg = 'black',
                               command=lambda: self.btn_return(master))
        return_btn.pack(side = tk.LEFT, padx = 10)
        self.pause_btn = tk.Button(master=self.button_frm, text='Pause', font = 'Arial 10', fg = 'black',
                               command=lambda: self.btn_pause(master))
        self.pause_btn.pack(side = tk.LEFT, padx = 10)
        quit_btn = tk.Button(master=self.button_frm, text='Quit', font='Arial 10 bold', fg = 'black', bg = 'red',
                             command=lambda: self.btn_quit(master))
        quit_btn.pack(side = tk.RIGHT, padx = 10)

    def btn_return(self, master):   # Return to start page
        self.plots_frm.destroy()
        self.button_frm.destroy()
        global defaults
        defaults = False
        master.switch_frame(StartPage)

    def btn_pause(self, master):    # Pause measurements
        if not self.paused:
            self.pause_btn.configure(text='Resume')
            self.paused = True
        else:
            self.pause_btn.configure(text='Pause')
            self.paused = False

    def btn_quit(self, master):     # Quit the program
        app.destroy()
        sys.exit()

    # SIGNIFICANT FIGURES ROUNDING
    def sigfig(self, x, n):  # Rounds number x to n significant figures
        if type(n) != int:
            #print('Error: non-integer significant figures')
            return x
        else:
            sf1 = math.ceil(-np.log10(abs(x)))  # finds first sig fig
            y = round(x, sf1 + n - 1)
            return y

    # PHASE-SENSITIVE DEMODULATION
    def demodulate(self, time, V_A, V_B):
        global params

        # MOVING AVERAGE FUNCTION
        def moving_average(data_array, window_length):  # Takes a number (window_length) of values in each average to return an array
            moving_average_array = []       # Initialize an empty list to store moving averages
            i = 0
            while i < len(data_array) - window_length + 1:
                window_average = np.sum(data_array[i:i + window_length]) / window_length    # Calculate the average of current window
                moving_average_array.append(window_average) # Store the average of current window in moving average list
                i += 1                      # Shift window to right by one position
            return moving_average_array

        demod = np.multiply(V_A, V_B)   # Demodulate the signal
        samplespace = params.sampling_time/params.Nsamples # Time between samples
        avg_window = int(round(params.averaging_time/samplespace, 0))
        DC_avg = moving_average(demod, avg_window)  # Find the moving DC average

        # Cut values off both ends of the time array to give it the same length as DC_avg
        avgtime = time[int(math.floor((len(time) - len(DC_avg)) / 2)): int(math.floor((len(time) + len(DC_avg)) / 2))]
        # Set the y-limits for the demodulation plot
        avgymax = max((max(DC_avg), abs(min(DC_avg))))
        # DC average value over whole demodulated array
        DC_avg_val = np.ones((np.size(avgtime))) * np.average(DC_avg)

        return DC_avg, DC_avg_val, avgtime, avgymax

    # CREATE THE PLOTS
    def setup_plot(self, masterframe):
        ani_fig = plt.figure(figsize = (12,7))
        ani_axA = ani_fig.add_subplot(221) # Add a subplot to a 2-row, 2-column grid in position 1
        ani_axB = ani_fig.add_subplot(222) # Add a subplot to a 2-row, 2-column grid in position 2
        ani_axPSD = ani_fig.add_subplot(212) # Add a subplot to a 2-row, 1-column grid in position 2

        self.plots = []  # Empty array of plots
        self.plots.append(ani_axA)
        self.plots.append(ani_axB)
        self.plots.append(ani_axPSD)

        def animate(i):
            if not self.paused:
                # Pull data off scope (block mode)
                global params
                global scope
                channel_configs = [device.ChannelConfig('A', True, params.Amode, params.Vinput_rangeA),
                                   device.ChannelConfig('B', True, params.Bmode, params.Vinput_rangeB)]
                timebase_options = device.TimebaseOptions(max_time_interval=None,
                                                          no_of_samples=None,
                                                          min_collection_time=params.sampling_time,     # Time to sample for
                                                          oversample=1)     # Don't oversample

                # Create the data
                times, voltages, overflow_warnings = scope.capture_block(timebase_options, channel_configs)
                crop_factor = int(round(len(times)/params.Nsamples))
                time = times[::crop_factor]             # Crop it down to N data points
                if params.modulation_channel == 'A':    # Take account of 10-times probes
                    A_voltages = voltages['A'][::crop_factor]
                    B_voltages = voltages['B'][::crop_factor] * 10
                else:
                    A_voltages = voltages['A'][::crop_factor] * 10
                    B_voltages = voltages['B'][::crop_factor]
                PSD_voltages, PSD_voltage, PSD_time, PSDVmax = self.demodulate(time, A_voltages, B_voltages)
                self.data = [A_voltages,B_voltages,PSD_voltages]
                self.times = [time, time, PSD_time]

                global PSDVlim
                if i == 0:
                    PSDVlim = PSDVmax

                if PSDVlim < PSDVmax:   # If the max value from this sample is greater than any previous sample, change the limit
                    PSDVlim = PSDVmax

                formats = ['b-', 'r-', 'm-']
                for i, ax in enumerate(self.plots):  # Iterate any repeated characteristics
                    ax.clear()
                    ax.plot(self.times[i], self.data[i], formats[i])
                    ani_axPSD.plot(PSD_time, PSD_voltage, 'c--', label=f'DC Average: {self.sigfig(PSD_voltage[0],4)}')
                    ani_axPSD.axhline(y=0, color='k', linewidth=1)
                    ax.set_xlabel('Time (s)')
                    ax.set_ylabel('Voltage (V)')
                    ax.set_xlim(min(time), max(time))

                ani_axA.set_ylim(-params.Vinput_rangeA, params.Vinput_rangeA)
                ani_axB.set_ylim(-params.Vinput_rangeB, params.Vinput_rangeB)

                ani_axPSD.set_title('Demodulated signal')
                ani_axPSD.set_xlim(min(PSD_time), max(PSD_time))
                ani_axPSD.set_ylim(-1.1*PSDVlim, 1.1*PSDVlim)
                ani_axPSD.legend(loc = 'upper center', framealpha = 1)

                # Check if any channels are going overrange
                if 'A' in overflow_warnings:
                    ani_axA.set_title('Channel A \nOVER RANGE')
                else:
                    ani_axA.set_title('Channel A')
                if 'B' in overflow_warnings:
                    ani_axB.set_title('Channel B \nOVER RANGE')
                else:
                    ani_axB.set_title('Channel B')

                ani_fig.tight_layout()
            else:
                sleep(0.2)

        canvas = FigureCanvasTkAgg(ani_fig, master=masterframe)
        canvas.draw()
        canvas.get_tk_widget().pack()
        canvas._tkcanvas.pack()
        app.ani = animation.FuncAnimation(ani_fig, animate, frames=None, repeat=False, cache_frame_data=False)

# PROTOCOL IF WINDOW IS CLOSED USING X IN CORNER
def on_quit():
    print('um')
    app.destroy()
    sys.exit()

# RUN THE SCRIPT
if __name__ == "__main__":
    global app
    app = gui()
    app.title("PicoScope PSD")
    app.protocol('WM_DELETE_WINDOW', on_quit)
    app.mainloop()
    sys.exit()
