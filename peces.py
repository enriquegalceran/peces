#####7segments
import RPi.GPIO as GPIO
import dht11
import time
import os
from Salida_limpia import mostrarresultados
from parametros import generate_parameter_list

# Define GPIO INPUT model
GPIO.setmode(GPIO.BCM)
GPIO.cleanup()

# temp and button pins
temp_pin = 19
button_pin = 26

GPIO.setup(temp_pin, GPIO.IN)
GPIO.setup(button_pin, GPIO.IN)

# segment pins
segments_pins_board = (11, 4, 16, 8, 7, 10, 18, 25)
for segment in segments_pins_board:
    GPIO.setup(segment, GPIO.OUT)
    GPIO.output(segment, 0)

# digit cathodes pins
digit_pins_board = (22, 27, 17, 12)
for digit in digit_pins_board:
    GPIO.setup(digit, GPIO.OUT)
    GPIO.output(digit, 1)

# Numbers and characters 7 segment layout dictionary
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

instance = dht11.DHT11(pin=temp_pin)
temp_readout_time = 2  # seconds
t_0 = time.time()
usb_directory = '/media/pi'



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


str_th, temp_hum = read_th((0, 0))
mode = 1
pressed = 0
button = 0

f = open("Data.txt", "w+")
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


finally:
    print('\n\nGoodbye!\n')
    GPIO.cleanup()
    f.close()
