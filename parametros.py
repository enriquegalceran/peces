"""
    System to define parameters using a file named param.txt
"""


def param_init():
    """
    Define initial parameters, which include default values which will be used

    :return:
    """
    param_key = ['temp_pin', 'button1_pin',
                 'segments_pins_board',
                 'digit_pins_board',
                 'temp_readout_time', 'write_wait_time', 'data_filename', 'usb_name']
    param_tuple = [0, 0, 1, 1, 0, 0, 0, 0]
    param_int = [1, 1, 1, 1, 1, 1, 0, 0]
    param_default = [19, 26,
                     (11, 4, 16, 8, 7, 10, 18, 25),
                     (22, 27, 17, 12),
                     2, 300, 'Data.txt', 'KINGSTON']
    return param_key, param_tuple, param_int, param_default


def read_file(paramfilename):
    """
    Read file and generate a list separating each line in the file

    :param paramfilename:
    :return:
    """
    param = []
    with open(paramfilename) as fp:
        line = fp.readline()
        param.append(line)
        while line:
            line = fp.readline()
            param.append(line)
    return param


def clean_list(param):
    """
    clean list, removing spaces and '\n'

    :param param:
    :return:
    """
    param_clean = []
    for line in param:
        text = line
        text = text.replace(' ', '').replace('\n', '')
        param_clean.append(text)
    return param_clean


def generate_list(param_clean, param_key, param_tuple, param_int):
    """
    Generates a list following the parameters set in previous lists. This allows the use of tuples, int, string or empty
    values.

    :param param_clean:
    :param param_key:
    :param param_tuple:
    :param param_int:
    :return:
    """
    # Initialize final param_value list
    param_value = []
    for i in range(len(param_key)):
        param_value.append('')

    for parameter in param_clean:
        tmp = parameter.split('=')
        if len(tmp) == 2:
            if tmp[0] in param_key and tmp[1] is not '':
                idx = param_key.index(tmp[0])
                if param_tuple[idx]:
                    value = tmp[1].replace('(', '').replace(')', '').split(',')
                    value = map(int, value)
                    param_value[idx] = tuple(value)
                else:
                    if param_int[idx]:
                        param_value[idx] = int(tmp[1])
                    else:
                        param_value[idx] = tmp[1]
    return param_value


def fill_default_values(param_value, param_default):
    """
    If some value has not been set, it will fill it with the default value

    :param param_value:
    :param param_default:
    :return:
    """
    for item in range(len(param_value)):
        if param_value[item] is '':
            param_value[item] = param_default[item]
    return param_value


def generate_parameter_list(filename='param.txt'):
    """
    Main function. Generates a list of parameters following the values set in the file 'filename'. Values that don't fit
    the desired values aren't considered and ignored. Empty values will be filled with a default values, located inside
    the param_init() function.

    :param filename:
    :return:
    """
    param_key, param_tuple, param_int, param_default = param_init()
    parameters = read_file(filename)
    parameters = clean_list(parameters)
    parameters = generate_list(parameters, param_key, param_tuple, param_int)
    parameters = fill_default_values(parameters, param_default)
    return parameters


if __name__ == "__main__":
    parameter_list = generate_parameter_list()
    print(parameter_list)

