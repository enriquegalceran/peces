#####7segments
import RPi.GPIO as GPIO
import dht11
import time
import os
from Salida_limpia import mostrarresultados
from parametros import generate_parameter_list
from datetime import datetime

# ---------------- Define Initial Parameters ---------------#
parameter_list = generate_parameter_list('param.txt')       #
temp_pin = parameter_list[0]                                # Pin in GPIO associated to temperature
button1_pin = parameter_list[1]                             # Pin in GPIO associated with the button
segments_pins_board = parameter_list[2]                     # Pins in GPIO for segments
digit_pins_board = parameter_list[3]                        # Pins in GPIO for digit cathodes
temp_readout_time = parameter_list[4]                       # Time in seconds to wait before reading temp and humidity
write_wait_time = parameter_list[5]                         # Time in seconds to wait before updating file
data_filename = parameter_list[6]                           # Name of the datafile where we will save the history
usb_name = parameter_list[7]                                # USB where data will be transferred
usb_directory = '/media/pi'                                 # Raspberry Pi's directory where the USB will appear
file_directory = '/home/pi'                                 # Raspberry Pi's directory where the file is
# ----------------------------------------------------------#

# set up GPIO
GPIO.setmode(GPIO.BCM)
GPIO.cleanup()
GPIO.setup(temp_pin, GPIO.IN)
GPIO.setup(button1_pin, GPIO.IN)

for segment in segments_pins_board:
    GPIO.setup(segment, GPIO.OUT)
    GPIO.output(segment, 0)

digit_pins_board = (22, 27, 17, 12)
for digit in digit_pins_board:
    GPIO.setup(digit, GPIO.OUT)
    GPIO.output(digit, 1)

# Numbers and characters 7 segment layout dictionary
instance = dht11.DHT11(pin=temp_pin)
segment_dictionary = {' ': (0, 0, 0, 0, 0, 0, 0),
                      '0': (1, 1, 1, 1, 1, 1, 0),
                      '1': (0, 1, 1, 0, 0, 0, 0),
                      '2': (1, 1, 0, 1, 1, 0, 1),
                      '3': (1, 1, 1, 1, 0, 0, 1),
                      '4': (0, 1, 1, 0, 0, 1, 1),
                      '5': (1, 0, 1, 1, 0, 1, 1),
                      '6': (1, 0, 1, 1, 1, 1, 1),
                      '7': (1, 1, 1, 0, 0, 0, 0),
                      '8': (1, 1, 1, 1, 1, 1, 1),
                      '9': (1, 1, 1, 1, 0, 1, 1),
                      'a': (1, 1, 1, 0, 1, 1, 1),
                      'b': (0, 0, 1, 1, 1, 1, 1),
                      'c': (1, 0, 0, 1, 1, 1, 0),
                      'd': (0, 1, 1, 1, 1, 0, 1),
                      'e': (1, 0, 0, 1, 1, 1, 1),
                      'f': (1, 0, 0, 0, 1, 1, 1),
                      'g': (1, 0, 1, 1, 1, 1, 0),
                      'h': (0, 1, 1, 0, 1, 1, 1),
                      'i': (0, 0, 0, 0, 1, 1, 0),
                      'j': (0, 1, 1, 1, 0, 0, 0),
                      'l': (0, 0, 0, 1, 1, 1, 0),
                      'n': (0, 0, 1, 0, 1, 0, 1),
                      'o': (0, 0, 1, 1, 1, 0, 1),
                      'p': (1, 1, 0, 0, 1, 1, 1),
                      'r': (0, 0, 0, 0, 1, 0, 1),
                      's': (1, 0, 1, 1, 0, 1, 1),
                      't': (0, 0, 0, 1, 1, 1, 1),
                      'u': (0, 1, 1, 1, 1, 1, 0),
                      'y': (0, 1, 1, 1, 0, 1, 1),
                      }


def read_th(temperature_humidity_old):
    readout = instance.read()
    if readout.is_valid():
        print('good read')
        text_output = str(str(readout.temperature) + str(readout.humidity)).rjust(4)
        return text_output, (readout.temperature, readout.humidity)
    else:
        print('Error in readout!')
        return str(str(temperature_humidity_old[0]) + str(temperature_humidity_old[1])).rjust(4), \
            temperature_humidity_old


def show_screen(text, dictionary, digit_pins=digit_pins_board, segments_pins=segments_pins_board):
    for dig in range(4):
        for segm in range(0, 7):
            GPIO.output(segments_pins[segm], dictionary[text[dig]][segm])
        if dig == 1:
            GPIO.output(segments_pins[-1], 1)
        else:
            GPIO.output(segments_pins[-1], 0)
        GPIO.output(digit_pins[dig], 0)
        time.sleep(0.005)
        GPIO.output(digit_pins[dig], 1)


def set_new_parameters(parameter_filename):
    global temp_pin, button1_pin, segments_pins_board, digit_pins_board, temp_readout_time, data_filename, \
        usb_directory, usb_name

    param_list = generate_parameter_list(parameter_filename)  #
    temp_pin = param_list[0]  # Pin in GPIO associated to temperature
    button1_pin = param_list[1]  # Pin in GPIO associated with the button
    segments_pins_board = param_list[2]  # Pins in GPIO for segments
    digit_pins_board = param_list[3]  # Pins in GPIO for digit cathodes
    temp_readout_time = param_list[4]  # Time in seconds to wait before reading temp and humidity
    data_filename = param_list[5]  # Name of the datafile where we will save the history
    usb_directory = '/media/pi'  # Raspberry Pi's directory where the USB will appear
    usb_name = 'KINGSTON'  # USB where data will be transferred


def add_line_to_file(filename, newline):
    with open(filename, 'a') as file:
        file.write(newline + '\n')


str_th, temp_hum = read_th((0, 0))
mode = 1
pressed = 0
button = 0
t_0 = time.time()
t_write_0 = -999
folder_directory = os.listdir(file_directory)
print(folder_directory)

# Check if the datafile is there.
if data_filename in folder_directory:
    add_line_to_file(data_filename, '\n\n--- New Execution ---\n')
else:
    open(data_filename, 'a').close()
    add_line_to_file(data_filename, '--- New Execution ---\n')


# Temporal para hacer la prueba
write_wait_time = 5
try:
    while True:
        t_1 = time.time()
        # Check if enough time has passed since last readout
        if t_1 - t_0 > temp_readout_time:
            t_0 = t_1
            str_th, temp_hum = read_th(temp_hum)

        # Change Mode
        button = GPIO.input(26)
        if pressed == 0 and button == 1:
            time.sleep(0.1)
            mode = (mode + 1) % 2
            pressed = 1
        elif button == 0:
            pressed = 0

        # Temperature/humidity mode
        if mode == 0:
            show_screen(str_th, segment_dictionary)

        # Time mode
        elif mode == 1:
            n = time.ctime()[11:13] + time.ctime()[14:16]
            str_time = str(n).rjust(4)
            show_screen(str_time, segment_dictionary)

        # write
        t_write_1 = time.time()
        # Check if enough time has passed since last readout
        if t_write_1 - t_write_0 > write_wait_time:
            now = datetime.now()
            # dd/mm/YY H:M:S
            dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
            print('Updating DataFile')
            t_write_0 = t_write_1
            add_line_to_file(data_filename, '{} | Temp= {} | Hum= {}'.format(dt_string, temp_hum[0], temp_hum[1]))



        if False: # Pulsado tecla de resetear parametros
            update_parameter = None



finally:
    print('\n\nGoodbye!\n')
    GPIO.cleanup()
