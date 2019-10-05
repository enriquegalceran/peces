# -*- coding: utf-8 -*-

"""
    Proyecto de Unidad de Limpieza y Tratamiento de Resultados Observacionales Nativo (Proyecto U.L.T.R.O.N.)

    Recoge los resultados de observaciones del instrumento CAFOS y los reduce correctamente.
    Desarrollado por Enrique Galceran García
"""

from astropy.io import fits
from astropy.time import Time
import numpy as np
import pandas as pd
import numpy.ma as ma
import matplotlib.pyplot as plt
import os
import csv
import argparse
import time
import datetime
import warnings
import json

from .Salida_limpia import mostrarresultados, stdrobust
from .IMGPlot import imgdibujar, limites_imagen
from .str2datetime import str2datetime


def tuple2coordinates(tupla):
    return tupla[0], tupla[1], tupla[2], tupla[3]


def save_file_csv(csvfile, res):
    with open(csvfile, "w") as output:
        writer = csv.writer(output, lineterminator='\n')
        for val in res:
            writer.writerow([val])


def show_dictionary(dictionary_name):
    print("Showing dictionary: ")
    for x, y in dictionary_name.items():
        print(x, y)


def save_json(variable, filename):
    json_f = json.dumps(variable, indent=2, sort_keys=True)
    f = open(filename, "w")
    f.write(json_f)
    f.close()


def load_json(filename='Dic_filtro.json'):
    with open(filename) as json_file:
        data = json.load(json_file)
    return data


def obtain_files_lists(path_):
    file_list = []
    for file in os.listdir(path_):
        if file.endswith(".fits"):
            file_list.append(os.path.join(path_, file))
    return file_list


def read_dictionary(name, dictionary_, filename='Dic_filtro.json'):
    """
    Busca en el diccionario el número del filtro buscado. Si no existe, crea uno nuevo

    :param name:
    :param dictionary_:
    :param filename:
    :return:
    """
    if name in dictionary_:
        return dictionary_[name]
    else:
        len_dic = len(dictionary_)
        dictionary_[name] = len_dic
        print('New dictionry entry: Filter {0} - Index {1:03d}'.format(name, dictionary_[name]))
        save_json(dictionary_, filename)
        return dictionary_[name]


def obtain_coordinates_ccd(image_, mypath_=False):
    """
    Given a string following the CAFOS structure, obtains the coordinates of the ccd (x1, x2, y1 and y2) as a tuple.

    :param image_:
    :param mypath_:
    :return:
    """
    if mypath_:
        str_ccdsec = fits.open(mypath_ + image_)[0].header['CCDSEC']
        longitud = len(str_ccdsec)
    else:
        str_ccdsec = image_
        longitud = len(str_ccdsec)+1

    comma = None
    for x in range(1, longitud):
        if str_ccdsec[x] == ',':
            comma = x
            break
    if comma is None:
        raise ValueError('comma not defined!')

    colon = None
    for x in range(comma + 1, longitud):
        if str_ccdsec[x] == ':':
            colon = x
            break
    if colon is None:
        raise ValueError('colon not defined!')

    coma2 = None
    for x in range(colon + 1, longitud):
        if str_ccdsec[x] == ',':
            coma2 = x
            break
    if coma2 is None:
        raise ValueError('coma2 not defined!')

    x1 = int(str_ccdsec[1:comma])
    y1 = int(str_ccdsec[comma + 1: colon])
    x2 = int(str_ccdsec[colon + 1: coma2])
    y2 = int(str_ccdsec[coma2 + 1: -1])
    tupla_salida = (x1, x2, y1, y2)
    return tupla_salida


def read_list(archivo):
    with open(archivo, 'rt') as f:
        reader = csv.reader(f, delimiter=',')
        your_list = list(reader)
    your_list = [item for sublist in your_list for item in sublist]
    return your_list


def obtain_coordinates_2(lista, idx):
    x1 = lista[idx, 0]
    x2 = lista[idx, 1]
    y1 = lista[idx, 2]
    y2 = lista[idx, 3]
    output_tuple = (x1, x2, y1, y2)
    return output_tuple


def obtain_bias(dir_bias_, night, list_nights, list_bias, x1, x2, y1, y2, b1, b2,
                max_search=10, fill_bias=680, verbose=False):
    """
    Searches the list of available bias for a bias with the same night the photo was taken. Also checks size nd binning.
    If it finds such a bias, it will use that one directly. If it does not find one, it will look for a different bias
    which fits the desired parametres for other nights. If there aren't any, it generates a flat bias.

    :param dir_bias_:
    :param night:
    :param list_nights:
    :param list_bias:
    :param x1:
    :param x2:
    :param y1:
    :param y2:
    :param b1:
    :param b2:
    :param max_search:
    :param fill_bias:
    :param verbose:
    :return:
    """

    searched_bias_name = dir_bias_ + night +\
        "-{0:04d}_{1:04d}_{2:04d}_{3:04d}-B{4:02d}_{5:02d}.fits".format(x1, x2, y1, y2, b1, b2)
    exist = None
    take_another = False
    if searched_bias_name in list_bias:
        exist = True
    else:
        position = list_nights.index(night)
        for i in range(1, max_search):
            for mult in [-1, 1]:
                indx = i * mult
                position_new = position + indx
                if position_new >= len(list_nights):
                    # Reached end of the year. Doesn't search for next year. Can be fixed with a bigger folder and
                    # placing every night in a same folder
                    break
                night = list_nights[position_new]
                searched_bias_new = dir_bias_ + night + \
                    "-{0:04d}_{1:04d}_{2:04d}_{3:04d}-B{4:02d}_{5:02d}.fits"\
                    .format(x1, x2, y1, y2, b1, b2)

                if searched_bias_new in list_bias:
                    searched_bias_name = searched_bias_new
                    exist = True
                    take_another = True

    if exist:
        if verbose:
            if take_another:
                print('Bias taken from another night. The taken bias is:')
                print(searched_bias_name)
            else:
                print('Bias exists.')
        searched_bias = fits.getdata(searched_bias_name, ext=0)
    else:
        if verbose:
            print('There are no nearby bias available. A filled bias has been generated.')
        naxis1_expected, naxis2_expected = obtain_naxis((x1, x2, y1, y2), b2, b1)

        searched_bias = np.full((naxis1_expected, naxis2_expected), fill_bias, dtype=float)

    return exist, searched_bias_name, searched_bias


def obtain_flats(dir_flats_, night_, list_nights, list_flats, x1, x2, y1, y2, b1, b2, id_filter_, max_search=10,
                 fill_flat=1, verbose=False):
    """
    Simillarly to obtain_bias, searches the list of available flats for a flat with the same night the photo was taken.
    Also checks size nd binning. If it finds such a bias, it will use that one directly. If it does not find one, it
    will look for a different bias which fits the desired parametres for other nights. If there aren't any, it generates
    a flat bias.

    :param dir_flats_:
    :param night_:
    :param list_nights:
    :param list_flats:
    :param x1:
    :param x2:
    :param y1:
    :param y2:
    :param b1:
    :param b2:
    :param id_filter_:
    :param max_search:
    :param fill_flat:
    :param verbose:
    :return:
    """

    searched_flat_name = dir_flats_ + night_ + "-{0:04d}_{1:04d}_{2:04d}_{3:04d}-B{4:02d}_{5:02d}-F{6:03d}.fits" \
        .format(x1, x2, y1, y2, b1, b2, id_filter_)
    exist = False
    take_another = False

    if searched_flat_name in list_flats:
        exist = True
    else:
        position = list_nights.index(night_)
        for i in range(1, max_search):
            if exist:
                break
            for mult in [-1, 1]:
                indx = i * mult
                position_new = position + indx
                if position_new >= len(list_nights):
                    # Reached end of the year. Doesn't search for next year. Can be fixed with a bigger folder and
                    # placing every night in a same folder
                    break
                night_ = list_nights[position_new]
                searched_flat_new = dir_flats_ + night_ + \
                    "-{0:04d}_{1:04d}_{2:04d}_{3:04d}-B{4:02d}_{5:02d}-F{6:03d}.fits" \
                    .format(x1, x2, y1, y2, b1, b2, id_filter_)

                if searched_flat_new in list_flats:
                    searched_flat_name = searched_flat_new
                    exist = True
                    take_another = True
                    break

    if exist:
        if verbose:
            if take_another:
                print('Flat taken from another night. The taken flat is:')
                print(searched_flat_name)
            else:
                print('Flat exists.')
        searched_flat = fits.getdata(searched_flat_name, ext=0)
    else:
        if verbose:
            print('There are no nearby flats available. A filled bias has been generated.')

        # searched_flat = None
        # if searched_flat is None:
        #     raise ValueError('No hay definido un valor por defecto para el flat')

        naxis1_expected, naxis2_expected = obtain_naxis((x1, x2, y1, y2), b2, b1)

        searched_flat = np.full((naxis1_expected, naxis2_expected), fill_flat, dtype=float)

    return exist, searched_flat_name, searched_flat


def create_circular_mask(h, w, center=None, radius=None):
    """
    Generates a circular mask based on height, width, centre and radius of an image.

    :param h:
    :param w:
    :param center:
    :param radius:
    :return:
    """

    if center is None:  # use the middle of the image
        center = [int(w/2), int(h/2)]
    if radius is None:  # use the smallest distance between the center and image walls
        radius = min(center[0], center[1], w-center[0], h-center[1])

    yc, xc = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((xc - center[0])**2 + (yc-center[1])**2)

    mask = dist_from_center <= radius
    return mask


def create_unique_list(dir_data, night, stuff_list, header, binning=False, filter_name=False):

    """
    Given a list of parameter, generates a series of lists:
        - List of files
        - List of unique values
        - For each value, how many times it appears
        - For every element of the list, generates a value of which 'unique value' appears
        - How many unique values does it have
        - Name of filter, if applicable

    :param dir_data:
    :param night:
    :param stuff_list:
    :param header:
    :param binning:
    :param filter_name:
    :return:
    """
    lista = []

    for image in stuff_list:
        lista.append(fits.open(dir_data + night + '/' + image)[0].header[header])
    unique_list, count_list = np.unique(lista, return_counts=True)  # Counting each list

    bin_sections = -1 * np.ones((len(unique_list), 2), dtype=int)
    indx_stuff = np.zeros(len(stuff_list), dtype=int)
    name_filter = []

    for i in range(len(stuff_list)):
        for j in range(len(unique_list)):
            if lista[i] == unique_list[j]:
                indx_stuff[i] = j
                if binning:
                    bin_sections[j, 0] = int(fits.open(dir_data + night + '/' + stuff_list[i])[0].header['ccdbinX'])
                    bin_sections[j, 1] = int(fits.open(dir_data + night + '/' + stuff_list[i])[0].header['ccdbinY'])
                if filter_name:
                    name_filter.append(fits.open(dir_data + night + '/' + stuff_list[i])[0].header['INSFLNAM'])
                break

    return lista, unique_list, count_list, indx_stuff, bin_sections, name_filter


def most_probable_image(archivo):
    image_data = fits.getdata(archivo, ext=0)
    v_median = np.median(image_data)
    v_sd = stdrobust(image_data)
    if v_sd < 50 and v_median < 900:
        probable = 0  # Bias
    elif v_sd < 500:
        probable = 1  # Arc
    else:
        probable = 2  # Flat

    return probable


def checking(file_, descriptor, descriptor2=None, descriptor3=None, verbose=False):
    if any([descriptor2, descriptor3]):
        descriptortemp = descriptor + descriptor2 + descriptor3
        coincides = True
        for texto in descriptortemp:
            if texto in fits.open(file_)[0].header['OBJECT']:
                # Some coincide without being science images
                coincides = False
                break
    else:
        coincides = False
        for texto in descriptor:
            if texto in fits.open(file_)[0].header['OBJECT']:
                # Both are the same
                coincides = True
                break
    if verbose:
        print(file_, fits.open(file_)[0].header['OBJECT'], fits.open(file_)[0].header['imagetyp'], coincides)

    # Maybe add something that checks for test

    return coincides


def file_list_2(path_, desc_bias, desc_flats, desc_arc, verbose=False, calysci=True):
    lista_cal = []
    lista_sci = []
    lista_misc = []
    lista_bias = []
    lista_flat = []
    lista_science = []
    lista_arc = []
    lista_else = []
    lista_files = []
    lista_wrong = []

    # Listing every file in the folder
    for file in os.listdir(path_):
        if file.endswith(".fits"):
            if calysci:
                if '-cal-' in file:
                    lista_cal.append(file)
                elif '-sci-' in file:
                    lista_sci.append(file)
                else:
                    lista_misc.append(file)
            lista_files.append(os.path.join(path_, file))

            # Splitting by IMAGETYP
            types = fits.open(path_ + file)[0].header['IMAGETYP'].strip()
            if not fits.open(path_ + file)[0].header['OBJECT'] == 'Test':  # Check it's not a test file
                if types == 'BIAS' or types == 'bias':
                    coincides = checking(path_ + file, desc_bias, verbose=verbose)
                    if coincides:
                        lista_bias.append(file)
                    else:
                        lista_wrong.append(file)

                elif types == 'flat':
                    coincides = checking(path_ + file, desc_flats, verbose=verbose)
                    if coincides:
                        lista_flat.append(file)
                    else:
                        lista_wrong.append(file)

                elif types == 'arc':
                    coincides = checking(path_ + file, desc_arc, verbose=verbose)
                    if coincides:
                        lista_arc.append(file)
                    else:
                        lista_wrong.append(file)

                elif types == 'science':
                    coincides = checking(path_ + file, desc_bias, desc_flats, desc_arc, verbose=verbose)
                    if coincides:
                        lista_science.append(file)
                    else:
                        lista_wrong.append(file)

                else:
                    lista_else.append(file)

    for file in lista_wrong:

        if calysci:
            if file in lista_sci:
                lista_science.append(file)

        probable = most_probable_image(path_ + file)

        if probable == 0:
            lista_bias.append(file)
        elif probable == 1:
            lista_arc.append(file)
        elif probable == 2:
            lista_flat.append(file)

    return lista_bias, lista_flat, lista_arc, lista_science, lista_files, lista_wrong


def create_list_cal_and_sci(lista_nights_, dir_lists_, dir_data_, desc_bias, desc_flats, desc_arc, verbose, calysci):
    i = 0
    for night in lista_nights_:
        i += 1
        if night not in os.listdir(dir_lists_):
            os.mkdir(dir_lists_ + night + '/')

            path_ = dir_data_ + night + '/'
            l_bias, l_flat, l_arc, l_ciencia, l_archivos, l_falla = file_list_2(path_,
                                                                                desc_bias,
                                                                                desc_flats,
                                                                                desc_arc,
                                                                                verbose,
                                                                                calysci)

            mostrarresultados(['Bias', 'Flat', 'Arc', 'Ciencia', 'Falla'],
                              [len(l_bias), len(l_flat), len(l_arc), len(l_ciencia), len(l_falla)],
                              titulo=night, contador=i, valor_max=len(lista_nights_))

            # save_file_csv(dir_lista_ + night + '/' + 'CAL.csv', cal)
            save_file_csv(dir_lists_ + night + '/' + 'SCI.csv', l_ciencia)
            save_file_csv(dir_lists_ + night + '/' + 'ARC.csv', l_archivos)
            save_file_csv(dir_lists_ + night + '/' + 'LBias.csv', l_bias)
            save_file_csv(dir_lists_ + night + '/' + 'LFlat.csv', l_flat)
            save_file_csv(dir_lists_ + night + '/' + 'LArc.csv', l_arc)
            save_file_csv(dir_lists_ + night + '/' + 'Errores.csv', l_falla)


def add_to_file(file, what):
    open(file, 'a+').write(what)


def slicing_data(slicing_push, size_da, size_mb):
    if slicing_push == "NW":
        s1 = size_da[0] - size_mb[0]
        s2 = size_da[0]
        s3 = 0
        s4 = size_mb[1]
    elif slicing_push == "SE":
        s1 = 0
        s2 = size_mb[0]
        s3 = size_da[1] - size_mb[1]
        s4 = size_da[1]
    elif slicing_push == "NE":
        s1 = size_da[0] - size_mb[0]
        s2 = size_da[0]
        s3 = 0
        s4 = size_mb[1]
    else:
        s1 = 0
        s2 = size_mb[0]
        s3 = 0
        s4 = size_mb[1]
    return s1, s2, s3, s4


def obt_naxis(x2, x1, binn):
    naxis = int((x2 - x1 + 1) / binn)
    if (x2 - x1 + 1) % binn != 0:
        naxis += 1
    return naxis


def obtain_naxis(coordenadas_tupla, binn_1, binn_2=None):
    if binn_2 is None:
        binn_2 = binn_1

    x1, x2, y1, y2 = tuple2coordinates(coordenadas_tupla)

    naxis1 = obt_naxis(y2, y1, binn_1)
    naxis2 = obt_naxis(x2, x1, binn_2)

    return naxis1, naxis2
########################################################################################################################


def join_bias_images(night, unique_sections_, coordinates_sections_, sections_count_, indx_section_,
                     bin_sections_, dir_bias_, dir_data_, list_bias_,
                     interactive=False, cut_extra=False, verbose_imagen=False):
    list_element = None
    for section in range(len(unique_sections_)):
        print('section: ' + str(section))
        coordinates_image = obtain_coordinates_2(coordinates_sections_, section)
        x1, x2, y1, y2 = tuple2coordinates(coordinates_image)

        # Obtain Binning
        ccdbinx = bin_sections_[section, 1]
        ccdbiny = bin_sections_[section, 0]

        naxis1_expected, naxis2_expected = obtain_naxis(coordinates_image, ccdbinx, ccdbiny)

        master_biases = np.zeros((sections_count_[section],
                                  naxis1_expected,
                                  naxis2_expected), dtype=float)

        index0 = 0
        slicing_push = False
        header = None
        for image in range(len(list_bias_)):
            if indx_section_[image] == section:
                image_file = dir_data_ + night + '/' + list_bias_[image]
                image_data = fits.getdata(image_file, ext=0)
                if image_data[:, :].shape == master_biases[index0, :, :].shape:
                    master_biases[index0, :, :] = image_data[:, :]                     # Juntar
                    if index0 == 0:
                        header = fits.open(image_file)[0].header
                else:
                    size_mb = master_biases[index0, :, :].shape
                    size_da = image_data[:, :].shape
                    if cut_extra:
                        if not slicing_push:
                            warnings.warn("Sizes Incompatible!")
                            print("Sizes incompatible:")
                            print("Data size: " + str(size_da))
                            print("Master Bias size: " + str(size_mb) + "\n")
                            slicing_push = (input("Slicing fits. "
                                                  "Select side towards to push (SW), SE, NE, NW: ") or "SW")
                        s1, s2, s3, s4 = slicing_data(slicing_push, size_da, size_mb)
                        master_biases[index0, :, :] = image_data[s1:s2, s3:s4]
                    else:
                        warnings.warn("Sizes Incompatible!")
                        print("Sizes incompatible:")
                        print("Data size: " + str(size_da))
                        print("Master Bias size: " + str(size_mb) + "\n")
                        print("Skipping current Master Bias.")
                        print("Consider using slicing with '--recortar'. ")
                        input("Press Enter to continue...")
                        break
                index0 += 1

        master_bias_colapsed = np.median(master_biases, axis=0)
        if verbose_imagen:
            imgdibujar(master_bias_colapsed, verbose_=1)
            plt.show()
        file_name = night + "-{0:04d}_{1:04d}_{2:04d}_{3:04d}-B{4:02d}_{5:02d}.fits".format(x1, x2, y1, y2,
                                                                                            ccdbinx, ccdbiny)

        mostrarresultados(['N', 'ccdbinY', 'ccbinX', 'A', 'B', '-1'],
                          [len(indx_section_[indx_section_ == section]), ccdbiny, ccdbinx,
                           naxis1_expected, naxis2_expected, file_name],
                          titulo='Bias Realizado')

        if header is None:
            raise ValueError('No se ha creado correctamente la cabecera')

        masterbias_header = header.copy()

        ahora = datetime.datetime.now()
        ahora_dt = Time(ahora, format='datetime', scale='utc')
        masterbias_header.add_history('Realizado Bias a partir de ' +
                                      str(len(indx_section_[indx_section_ == section])) + ' imagenes. | ' +
                                      str(ahora_dt)[:19])
        numero_bias = 0
        for bias_raw in range(len(list_bias_)):
            if indx_section_[bias_raw] == section:
                numero_bias += 1
                masterbias_header.add_history(str(numero_bias) + ': ' + list_bias_[bias_raw])
                print(str(numero_bias) + ': ' + list_bias_[bias_raw])

        masterbias_final = fits.PrimaryHDU(master_bias_colapsed.astype(np.float32), masterbias_header)

        masterbias_final.writeto(dir_bias_ + file_name, overwrite=True)

        if verbose_imagen:
            coord_lim = limites_imagen(*coordinates_image)
            imgdibujar(master_bias_colapsed, *coordinates_image, *coord_lim, verbose_=1)

        if interactive:
            input("Press Enter to continue...")

        # Add to table
        isot_time = masterbias_header['DATE']
        time_ = Time(isot_time, format='isot', scale='utc')
        if list_element is None:
            list_element = pd.DataFrame([[naxis1_expected, naxis2_expected,
                                          x1, x2, y1, y2,
                                          ccdbinx, ccdbiny,
                                          file_name, night,
                                          time_, time_.jd]],
                                        columns=['Naxis1', 'Naxis2',
                                                 'x1', 'x2', 'y1', 'y2',
                                                 'Binning1', 'Binning2',
                                                 'nombre_archivo', 'noche',
                                                 'tiempo_astropy', 'julian'])
        else:
            list_element_ = pd.DataFrame([[naxis1_expected, naxis2_expected,
                                           x1, x2, y1, y2,
                                           ccdbinx, ccdbiny,
                                           file_name, night,
                                           time_, time_.jd]],
                                         columns=['Naxis1', 'Naxis2',
                                                  'x1', 'x2', 'y1', 'y2',
                                                  'Binning1', 'Binning2',
                                                  'nombre_archivo', 'noche',
                                                  'tiempo_astropy', 'julian'])
            list_element = pd.concat([list_element, list_element_], ignore_index=True)

    return list_element


def make_master_bias(list_nights, dir_lists, dir_data, dir_bias, interactive, cut_reshape,
                     verbose_imagen=False):
    i_noche = 0
    df_bias = None
    for night in list_nights:
        i_noche += 1
        print('=== NIGHT ' + night + ' - (' + str(i_noche) + '/' + str(len(list_nights)) + ') ===')

        list_bias = read_list(dir_lists + night + '/' + 'LBias.csv')

        sections, unique_section, section_counting, section_index, bin_sections, _ = create_unique_list(
            dir_data, night, list_bias, header='CCDSEC', binning=True
        )

        sections_coordinates = np.zeros((len(unique_section), 4), dtype=int)
        for i in range(len(unique_section)):
            unique_coordinates = obtain_coordinates_ccd(unique_section[i])
            sections_coordinates[i, :] = [*unique_coordinates]

        df_bias_ = join_bias_images(night, unique_section, sections_coordinates, section_counting,
                                    section_index, bin_sections, dir_bias, dir_data, list_bias,
                                    interactive=interactive, cut_extra=cut_reshape,
                                    verbose_imagen=verbose_imagen)
        if df_bias is None:
            df_bias = df_bias_
        else:
            df_bias = pd.concat([df_bias, df_bias_], ignore_index=True)

    return df_bias


def join_flat_images(nights, unique_sections_, sections_coordinates_, indx_section_,
                     dir_bias_, dir_data_, dir_flats_, list_flats_, list_nights_, lista_bias,
                     verbose=0, interactive=False, verbose_images=False):
    """
    Combine multiple flats to create a master flat.

    :param nights:
    :param unique_sections_:
    :param sections_coordinates_:
    :param indx_section_:
    :param dir_bias_:
    :param dir_data_:
    :param dir_flats_:
    :param list_flats_:
    :param list_nights_:
    :param lista_bias:
    :param verbose:
    :param interactive:
    :param verbose_images:
    :return:
    """

    # Load Dictionary File
    dic_filter = load_json()
    element_list = None

    # Splitting by night
    for section in range(len(unique_sections_)):
        print('section: ' + str(section))
        coordinates_image = obtain_coordinates_2(sections_coordinates_, section)
        x1, x2, y1, y2 = tuple2coordinates(coordinates_image)

        # Create a list with every one that fits the index
        list_coincide = []
        for image in range(len(list_flats_)):
            if indx_section_[image] == section:
                list_coincide.append(list_flats_[image])

        filters, unique_filters, filter_counts, indx_filter, _, name_filter = create_unique_list(
            dir_data_, nights, list_coincide, header='INSFLID', binning=False, filter_name=True
        )

        # crete unique_names
        unique_names = []
        for i in range(len(unique_filters)):
            for j in range(len(filters)):
                if unique_filters[i] == filters[j]:
                    unique_names.append(name_filter[j])
                    break

        # For each element of unique_names, we search it's intex in the dictionary
        name_dictionary = []
        for i in unique_names:
            name_dictionary.append(read_dictionary(i, dic_filter))

        # Only necessary if the ISN'T a grisma, so we search what it's grisma is being used
        grisma_number = int(fits.open(dir_data_ + nights + '/' + list_coincide[0])[0]
                            .header['insgrid'].replace(' ', '0')[6:8])

        if grisma_number == 11:  # If grisma_number==11, there is no grisma being used
            free_grisma = True
        else:
            free_grisma = False

        # We want to know what is the number of the filter used
        filter_number = np.zeros(len(unique_filters), dtype=int)
        p = 0
        for i in unique_filters:
            filter_number[p] = int(i[-2:].strip())
            p += 1

        # Separated by sections, so we separate by filter
        for filtro in range(len(unique_filters)):
            lista_actual = []
            binning2 = []
            binning1 = []

            # Binning is constant, so we split by different binning
            for i in range(len(list_coincide)):
                if indx_filter[i] == filtro:
                    lista_actual.append(list_coincide[i])
                    binning2.append(int(fits.open(dir_data_ + nights + '/' + list_coincide[i])[0].header['ccdbinX']))
                    binning1.append(int(fits.open(dir_data_ + nights + '/' + list_coincide[i])[0].header['ccdbinY']))

            bin1_unique, bin1_count = np.unique(binning1, return_counts=True)

            binning_always_coincide = True
            for i in range(len(lista_actual)):
                if not binning2 == binning1:
                    binning_always_coincide = False
                    break

            # If binning coincide
            master_flats_colapsed = None

            if binning_always_coincide:
                for binning in bin1_unique:
                    lista_actual_bin = []

                    # We have a list with individual binnings
                    for i in range(len(lista_actual)):
                        if binning1[i] == binning:
                            lista_actual_bin.append(lista_actual[i])

                    # Only has one bin and it matches, it's axis (image isn't stretched),
                    # we don't worry about it and generate the image
                    ccdbiny = binning
                    ccdbinx = binning

                    naxis1_expected, naxis2_expected = obtain_naxis(coordinates_image, ccdbinx, ccdbiny)

                    master_flats = np.zeros((filter_counts[filtro], naxis1_expected, naxis2_expected), dtype=float)

                    # Generate mask
                    real_center = [1075, 1040]
                    center = [real_center[0] - x1, real_center[1] - y1]
                    radius = 809  # con 810 tiene un pixel de borde
                    mask = create_circular_mask(naxis1_expected, naxis2_expected, center=center, radius=radius)

                    mostrarresultados(['N', 'ccdbiny', 'ccdbinx', 'A', 'B'],
                                      [len(lista_actual_bin), ccdbiny, ccdbinx,
                                       naxis1_expected, naxis2_expected],
                                      titulo='Filtro ' + str(name_dictionary[filtro]))

                    # Read for each element their values and we concatenate
                    index0 = 0
                    header_ = fits.open(dir_data_ + nights + '/' + lista_actual_bin[0])[0].header
                    for image in range(len(lista_actual_bin)):
                        image_file = dir_data_ + nights + '/' + lista_actual_bin[image]
                        image_data = fits.getdata(image_file, ext=0)

                        if image_data[:, :].shape == master_flats[index0, :, :].shape:
                            master_flats[index0, :, :] = image_data[:, :]                  # Juntar
                        else:
                            print('There is a problem with the image size')
                            input('Pause. Press Enter to to continue...')
                        index0 += 1

                    exist, searched_bias_name, searched_bias = obtain_bias(
                        dir_bias_, nights, list_nights_, lista_bias, x1, x2, y1, y2, ccdbinx, ccdbiny
                    )

                    for i in range(master_flats.shape[0]):
                        master_flats[i, :, :] = master_flats[i, :, :] - searched_bias

                    # Normalize
                    median_value = np.zeros(master_flats.shape[0], dtype=float)
                    for i in range(master_flats.shape[0]):

                        # If it has a grisma, we need to use the mask
                        if free_grisma:
                            x = ma.masked_array(master_flats[i, :, :], ~mask)
                            median_value[i] = ma.median(x)
                        else:
                            median_value[i] = np.median(master_flats[i, :, :])

                        if median_value[i] <= 0:
                            if median_value[i] == 0:
                                median_value[i] = 1
                            else:
                                median_value[i] = abs(median_value[i])

                        # Divide
                        master_flats[i, :, :] = np.true_divide(master_flats[i, :, :], median_value[i])

                    mean_value2 = np.zeros(master_flats.shape[0], dtype=float)
                    for i in range(master_flats.shape[0]):
                        mean_value2[i] = np.mean(master_flats[i, :, :], dtype=float)

                    # Collapse
                    master_flats_colapsed = np.median(master_flats, axis=0)

                    # Set corners as 1 to avoid dividing by 0
                    if free_grisma:
                        master_flats_colapsed[~mask] = 1

                    if verbose_images:
                        imgdibujar(master_flats_colapsed)
                        plt.show()

                    # Generate the name of the future flat
                    flat_file_name = nights +\
                        "-{0:04d}_{1:04d}_{2:04d}_{3:04d}-B{4:02d}_{5:02d}-F{6:03d}.fits"\
                        .format(x1, x2, y1, y2, ccdbinx, ccdbiny, int(name_dictionary[filtro]))

                    masterflats_header = header_.copy()
                    if masterflats_header['BLANK']:
                        del masterflats_header['BLANK']

                    now_time = datetime.datetime.now()
                    now_datetime = Time(now_time, format='datetime', scale='utc')
                    masterflats_header.add_history('Realizado Flat a partir de ' +
                                                   str(len(indx_section_[indx_section_ == section])) +
                                                   ' imagenes. | ' + str(now_datetime)[:19])
                    flat_number = 0
                    for flat_raw in range(len(list_flats_)):
                        if indx_section_[flat_raw] == section:
                            flat_number += 1
                            masterflats_header.add_history(str(flat_number) + ': ' + list_flats_[flat_raw])
                            print(str(flat_number) + ': ' + list_flats_[flat_raw])

                    if master_flats_colapsed is None:
                        raise ValueError('Collapsed flat created wrong.')
                    masterflats_final = fits.PrimaryHDU(master_flats_colapsed.astype(np.float32), masterflats_header)

                    masterflats_final.writeto(dir_flats_ + flat_file_name, overwrite=True)

                    # Add to table
                    isot_time = masterflats_header['DATE']
                    time_dataframe = Time(isot_time, format='isot', scale='utc')
                    if element_list is None:
                        element_list = pd.DataFrame([[naxis1_expected, naxis2_expected,
                                                      x1, x2, y1, y2,
                                                      ccdbinx, ccdbiny,
                                                      int(name_dictionary[filtro]),
                                                      free_grisma, grisma_number,
                                                      flat_file_name, nights,
                                                      time_dataframe, time_dataframe.jd]],
                                                    columns=['Naxis1', 'Naxis2',
                                                             'x1', 'x2', 'y1', 'y2',
                                                             'Binning1', 'Binning2',
                                                             'filtro',
                                                             'free_grisma', 'num_grisma',
                                                             'flat_file_name', 'noche',
                                                             'tiempo_astropy', 'julian'])
                    else:
                        elemento_lista_ = pd.DataFrame([[naxis1_expected, naxis2_expected,
                                                         x1, x2, y1, y2,
                                                         ccdbinx, ccdbiny,
                                                         int(name_dictionary[filtro]),
                                                         free_grisma, grisma_number,
                                                         flat_file_name, nights,
                                                         time_dataframe, time_dataframe.jd]],
                                                       columns=['Naxis1', 'Naxis2',
                                                                'x1', 'x2', 'y1', 'y2',
                                                                'Binning1', 'Binning2',
                                                                'filtro',
                                                                'free_grisma', 'num_grisma',
                                                                'flat_file_name', 'noche',
                                                                'tiempo_astropy', 'julian'])
                        element_list = pd.concat([element_list, elemento_lista_], ignore_index=True)

            else:
                raise ValueError("Binning don't fit on both axis")

            if verbose >= 1:
                coord_lim = limites_imagen(*coordinates_image)
                imgdibujar(master_flats_colapsed, *coordinates_image, *coord_lim, verbose_=1)

            if interactive:
                input("Press Enter to continue...")

    return element_list


def make_master_flat(list_nights, list_bias, dir_lists, dir_data, dir_bias, dir_flats,
                     verbose, interactive, verbose_imagen=False):

    """
    Main function to make master flats.

    :param list_nights:
    :param list_bias:
    :param dir_lists:
    :param dir_data:
    :param dir_bias:
    :param dir_flats:
    :param verbose:
    :param interactive:
    :param verbose_imagen:
    :return:
    """

    i_night = 0
    df_flat = None
    for night in list_nights[:]:
        i_night += 1
        print('=== NIGHT ' + night + ' - (' + str(i_night) + '/' + str(len(list_nights)) + ') ===')

        lista_flats = read_list(dir_lists + night + '/' + 'LFlat.csv')

        sections, unique_section, section_counting, section_index, bin_sections, _ = create_unique_list(
            dir_data, night, lista_flats, header='CCDSEC', binning=True
        )

        sections_coordinates = np.zeros((len(unique_section), 4), dtype=int)
        for i in range(len(unique_section)):
            unique_coordinates = obtain_coordinates_ccd(unique_section[i])
            sections_coordinates[i, :] = [*unique_coordinates]

        df_flat_ = join_flat_images(night, unique_section, sections_coordinates, section_index,
                                    dir_bias, dir_data, dir_flats, lista_flats, list_nights, list_bias,
                                    verbose=verbose, interactive=interactive, verbose_images=verbose_imagen)

        if df_flat is None:
            df_flat = df_flat_
        else:
            df_flat = pd.concat([df_flat, df_flat_], ignore_index=True)

    return df_flat


def reducing_images(list_nights, dir_listas, dir_datos, dir_bias, dir_flats, dir_reducc,
                    df_bias, df_flat, verbose=2, verbose_imagen=False):
    # Load lists
    element_list = None
    dic_filter = load_json()
    doesnt_exist = []
    total_science_images = 0
    saved_images = 0

    for night in list_nights:
        images_reduced_night = 0
        print(night)
        if night not in os.listdir(dir_reducc):
            os.mkdir(dir_reducc + night + '/')

        list_science = read_list(dir_listas + night + '/' + 'SCI.csv')
        doesnt_exist = (0, 0, len(list_science))

        secc, secc_unicas, secc_count, indice_secc, bin_secc, nombres_filtros = create_unique_list(
            dir_datos, night, list_science, header='CCDSEC', binning=True, filter_name=True
        )
        i_image = 0
        for image in range(len(list_science)):
            i_image += 1
            # Image name
            name_science = dir_datos + night + '/' + list_science[image]
            header = fits.open(name_science)[0].header

            # Coordinates and binning
            coordinates = obtain_coordinates_ccd(secc_unicas[indice_secc[image]])
            x1, x2, y1, y2 = tuple2coordinates(coordinates)
            binning = bin_secc[indice_secc[image]]
            naxis1_r = header['Naxis1']
            naxis2_r = header['Naxis2']

            naxis1_science, naxis2_science = obtain_naxis((x1, x2, y1, y2), binning[1], binning[0])

            # Check for overscan
            biassec = header['BIASSEC']
            coordinates_biassec = obtain_coordinates_ccd(biassec)
            naxis1_overscan, naxis2_overscan = obtain_naxis(coordinates_biassec, binning[1], binning[0])

            if coordinates_biassec[0] == 0 and coordinates_biassec[1] == 0:
                overscan = False
                x_1, x_2, y_1, y_2 = None, None, None, None
            else:
                print('There is overscan!')
                overscan = True

                # There is a need to cut the image. We generate the new section for the image
                if coordinates_biassec[0] > x2:
                    x_1 = 0
                    x_2 = naxis2_science
                    y_1 = 0
                    y_2 = naxis1_science
                elif coordinates_biassec[3] > y2:
                    x_1 = 0
                    x_2 = naxis2_science
                    y_1 = 0
                    y_2 = naxis1_science
                elif coordinates_biassec[1] < x1:
                    x_1 = naxis2_overscan + 1
                    x_2 = naxis2_science + naxis2_overscan
                    y_1 = 0
                    y_2 = naxis1_science
                elif coordinates_biassec[4] < y1:
                    x_1 = 0
                    x_2 = naxis2_science
                    y_1 = naxis1_overscan + 1
                    y_2 = naxis1_science + naxis2_overscan
                else:
                    raise ValueError('Where is the overscan?!')

            # We obtain the name of the filter and it's ID_filter
            filtername = nombres_filtros[image]
            id_filter = read_dictionary(filtername, dic_filter)

            # Julian Date
            isot_time = header['DATE']
            date_datetime = str2datetime(isot_time)
            time_isot = Time(isot_time, format='isot', scale='utc')

            if verbose > 0:
                mostrarresultados(['noche', 'naxis1_r', 'naxis2_r', 'naxis1_c', 'naxis2_c',
                                   'x1', 'x2', 'y1', 'y2',
                                   'binning', 'filtro', 'id_filtro',
                                   'FJM', 'Overscan',
                                   'imagen', 'nombre',
                                   'fecha', 'hora'],
                                  [night, naxis1_r, naxis2_r, naxis1_science, naxis2_science,
                                   x1, x2, y1, y2,
                                   binning, filtername, id_filter,
                                   time_isot.jd, overscan,
                                   image, list_science[image],
                                   date_datetime.date(), date_datetime.time()],
                                  contador=i_image, valor_max=len(list_science))

            # Buscamos el Bias
            # Seleccionamos los que coinciden con el tamaño teorico
            datat = df_bias[df_bias.Naxis1 == naxis1_science]
            datat = datat[datat.Naxis2 == naxis2_science]
            datat = datat[datat.Binning1 == binning[0]]
            datat = datat[datat.Binning2 == binning[1]]
            fechasjd_b = abs(datat.julian.values - time_isot.jd)
            if len(fechasjd_b) > 0:
                pos_min = np.argmin(fechasjd_b)
                nombre_bias_buscado = datat.iloc[pos_min]['nombre_archivo']
                if datat.iloc[pos_min]['noche'] != night:
                    if verbose > 1:
                        print('El bias mas cercano no es de esta noche.')

                bias_asociado = fits.getdata(dir_bias + nombre_bias_buscado, ext=0)
            else:
                nombre_bias_buscado = 'Ausente'
                relleno_b = 680
                if verbose > 1:
                    print('No se han encontrado bias de esa forma, se genera uno artificial.')
                bias_asociado = np.full((naxis1_science, naxis2_science), relleno_b, dtype=float)

            # Buscamos el Flat
            # Seleccionamos los que coinciden con el tamaño teorico
            dataf = df_flat[df_flat.Naxis1 == naxis1_science]
            dataf = dataf[dataf.Naxis2 == naxis2_science]
            dataf = dataf[dataf.Binning1 == binning[0]]
            dataf = dataf[dataf.Binning2 == binning[1]]
            dataf = dataf[dataf.filtro == id_filter]

            # Lista de los bias en función de la distancia en tiempo a la imagen de ciencia
            fechasjd_f = abs(dataf.julian.values - time_isot.jd)
            if len(fechasjd_f) > 0:
                pos_min = np.argmin(fechasjd_f)
                nombre_flat_buscado = dataf.iloc[pos_min]['nombre_archivo']
                if dataf.iloc[pos_min]['noche'] != night:
                    if verbose > 1:
                        print('Existe un Flat, pero no es de esta noche.')
                flat_asociado = fits.getdata(dir_flats + nombre_flat_buscado, ext=0)

            else:
                nombre_flat_buscado = 'Ausente'
                relleno_f = 1
                print('No se han encontrado flats, se genera uno artificialmente.')
                flat_asociado = np.full((naxis1_science, naxis2_science), relleno_f, dtype=float)

            # Obtenemos la información de la imagen de ciencia
            image_data = fits.getdata(name_science, ext=0)
            reducido_header = header.copy()

            # Comprobamos si hay o no overscan que haya que tener en cuenta
            if overscan:
                print('Recortamos los datos')
                image_data = image_data[y_1:y_2, x_1:x_2]

            ##############################################################################
            reducido_datos = (image_data - bias_asociado) / flat_asociado
            ##############################################################################

            if verbose_imagen:
                imgdibujar(reducido_datos)
                plt.show()

            if reducido_header['BLANK']:
                del reducido_header['BLANK']

            # Sacar si tiene o no un grisma
            numero_grisma = reducido_header['insgrid'].replace(' ', '0')[6:8]
            if numero_grisma == 11:  # Si numero_grisma==11, entonces no hay grisma de por medio
                free_grisma = True
            else:
                free_grisma = False

            if free_grisma:
                mask = create_circular_mask(naxis1_science, naxis2_science)  # Se puede cambiar el centro y radio
                reducido_datos[~mask] = 1

            ahora = datetime.datetime.now()
            ahora_dt = Time(ahora, format='datetime', scale='utc')
            reducido_header.add_history('Realizada la reduccion a partir de las siguientes imagenes: | ' +
                                        str(ahora_dt)[:19])
            reducido_header.add_history('Bias: ' + nombre_bias_buscado)
            reducido_header.add_history('Flat: ' + nombre_flat_buscado)

            # Guardamos la imagen
            reducido_final = fits.PrimaryHDU(reducido_datos.astype(np.float32), reducido_header)

            reducido_final.writeto(dir_reducc + night + '/' + 'r_' + list_science[image], overwrite=True)

            images_reduced_night += 1

            ahora = datetime.datetime.now()
            time_isot = Time(ahora, format='datetime', scale='utc')
            if element_list is None:
                element_list = pd.DataFrame([[naxis1_science, naxis2_science,
                                                x1, x2, y1, y2,
                                                binning[0], binning[1],
                                                id_filter,
                                                free_grisma, numero_grisma,
                                                'r_' + list_science[image],
                                                nombre_bias_buscado, nombre_flat_buscado,
                                                night, time_isot, time_isot.jd]],
                                              columns=['Naxis1', 'Naxis2',
                                                       'x1', 'x2', 'y1', 'y2',
                                                       'Binning1', 'Binning2',
                                                       'filtro',
                                                       'free_grisma', 'num_grisma',
                                                       'nombre_archivo',
                                                       'nombre bias', 'nombre flat',
                                                       'noche', 'Fecha realizacion', 'julian'])
            else:
                elemento_lista_ = pd.DataFrame([[naxis1_science, naxis2_science,
                                                 x1, x2, y1, y2,
                                                 binning[0], binning[1],
                                                 id_filter,
                                                 free_grisma, numero_grisma,
                                                 'r_' + list_science[image],
                                                 nombre_bias_buscado, nombre_flat_buscado,
                                                 night, time_isot, time_isot.jd]],
                                               columns=['Naxis1', 'Naxis2',
                                                        'x1', 'x2', 'y1', 'y2',
                                                        'Binning1', 'Binning2',
                                                        'filtro',
                                                        'free_grisma', 'num_grisma',
                                                        'nombre_archivo',
                                                        'nombre bias', 'nombre flat',
                                                        'noche', 'Fecha realizacion', 'julian'])
                element_list = pd.concat([element_list, elemento_lista_], ignore_index=True)

        # Al final de cada noche se hace el recuento
        doesnt_exist.append(doesnt_exist)
        saved_images += images_reduced_night
        print('verbosidad reduccion', verbose)
        if verbose > 0:
            mostrarresultados(['Imagenes reducidas', 'Acumulado'], [images_reduced_night, saved_images],
                              titulo="Reduccion de Imagenes",
                              contador=list_nights.index(night), valor_max=len(list_nights))

    mostrarresultados(list_nights, doesnt_exist)
    no_existen_2 = [0, 0]
    _ = element_list.to_csv('df_reducido.csv', index=None, header=True)
    for i in range(len(doesnt_exist)):
        total_science_images += doesnt_exist[i][2]
        no_existen_2[0] += doesnt_exist[i][0]
        no_existen_2[1] += doesnt_exist[i][1]

    print('Biases que no ha funcionado: ', no_existen_2[0], '| Flats que no ha funcionado: ', no_existen_2[1],
          '| Imagenes en total: ', total_science_images, '| Imagenes reducidas: ', saved_images)

    return saved_images


def pruebas_pandas(data):
    print(data.head())
    # print(data[data.noche == '170225_t2_CAFOS'])
    datat = data[data.Naxis1 == 2048]
    datat = datat[datat.Naxis2 == 1000]
    print(datat)
    print(data.nombre_archivo.values[0])


def decidir_repetir_calculos(norealizar, sirealizar, sujeto, dir_df, dir_sujetos):
    existe_bias = os.path.exists(dir_df + 'df_' + sujeto + '.csv')
    lista_existentes = os.listdir(dir_sujetos)
    if not any([norealizar, sirealizar]):
        if existe_bias:
            print('Se va a coger una version ya existente de los ' + sujeto)
            importado = pd.read_csv(dir_df + 'df_' + sujeto + '.csv')
            lista_teorica = importado.nombre_archivo.values.tolist()

            contador_t = 0
            contador_e = 0
            for i in range(len(lista_existentes)):
                if lista_teorica[i] in lista_existentes:
                    contador_t += 1
                if lista_existentes[i] in lista_teorica:
                    contador_e += 1

            if contador_e == len(lista_existentes) and contador_t == len(lista_teorica):
                print('Estan todos contabilizados, no hace falta volver a calcular los ' + sujeto)
                realizar = False
            else:
                print('No estan todos los ' + sujeto + ', se vuelven a calcular')
                realizar = True
        else:
            print('No existe el dataframe, se calcula el dataframe de los ' + sujeto)
            realizar = True
    elif norealizar:
        if existe_bias:
            print('Existe el dataframe de ' + sujeto + '. No se repite')
            realizar = False
        else:
            print('No existe el dataframe de ' + sujeto + ', se fuerza que se vuelva a hacer')
            realizar = True
    elif sirealizar:
        if existe_bias:
            print('Existe el dataframe de ' + sujeto + '. Se fuerza calcularlo de nuevo')
            realizar = True
        else:
            print('No existe el dataframe de ' + sujeto + ', Se fuerza que se haga')
            realizar = True
    else:
        raise ValueError('No cuadran las condiciones para los ' + sujeto + '!')

    return realizar


def main():

    # ---------------Valores por defecto-------------------------------------------
    default_dir_datos = '/media/enrique/TOSHIBA EXT/DataCAHA/CAFOS2017/'
    default_dir_bias = '/media/enrique/TOSHIBA EXT/CAHA/Biases/'
    default_dir_listas = '/media/enrique/TOSHIBA EXT/CAHA/Listas/'
    default_dir_flats = '/media/enrique/TOSHIBA EXT/CAHA/Flats/'
    default_dir_reduccion = '/media/enrique/TOSHIBA EXT/CAHA/Reduccion/'
    default_dir_dataframe = ''
    desc_bias = ['bias', 'Bias', 'BIAS']
    desc_flats = ['flats', 'FLATS', 'FLAT', 'Flats', 'Flat', 'flat', 'Skyflat', 'SDkyflat']
    desc_arc = ['arc', 'ARC']
    # -----------------------------------------------------------------------------

    parser = argparse.ArgumentParser(description="Bias and Flat calibration of CAFOS images")
    group = parser.add_mutually_exclusive_group()
    grupo2 = parser.add_mutually_exclusive_group()
    grupo3 = parser.add_mutually_exclusive_group()

    group.add_argument("-v", "--verbose", action="store_true")
    group.add_argument("-q", "--quiet", action="store_true")
    parser.add_argument("-db", "--dir_bias", default=default_dir_bias, type=str, help='Bias Directory')
    parser.add_argument("-df", "--dir_flats", default=default_dir_flats, type=str, help='Flats Directory')
    parser.add_argument("-dd", "--dir_datos", default=default_dir_datos, type=str, help='Data Directory')
    parser.add_argument("-dl", "--dir_listas", default=default_dir_listas, type=str, help='Lists Directory')
    parser.add_argument("-de", "--dir_reducc", default=default_dir_reduccion, type=str, help='Reducction Directory')
    parser.add_argument("-ddf", "--dir_dataf", default=default_dir_dataframe, type=str, help='DataFrame Directory')
    parser.add_argument("-vi", "--verboseimage", action="store_true", help="Mostrar Imagenes")
    parser.add_argument('--cmap', type=str, help="Colormap", default='hot')
    parser.add_argument("-i", "--interactive", action="store_true")
    parser.add_argument("--recortar", action="store_true", help="Activar el recorte de imagenes")
    parser.add_argument("--calysci", action="store_false",
                        help="Usar cuando los archivos no tienen '-cal-' y '-sci-'"
                             + "en el nombre para diferenciar entre calibración y ciencia.")
    parser.add_argument("-nr", "--noreducc", action="store_false", help="No realizar la reduccion.")

    grupo2.add_argument("-nb", "--nobias", action="store_true",
                        help="No realizar los Master Bias.")
    grupo2.add_argument("-sb", "--sibias", action="store_true",
                        help="Fuerza realizar los Master Bias.")

    grupo3.add_argument("-nf", "--noflat", action="store_true",
                        help="No realizar los Master Flat.")
    grupo3.add_argument("-sf", "--siflat", action="store_true",
                        help="Fuerza realizar los Master Flat.")

    args = parser.parse_args()

    if args.verbose:
        verbosidad = 2
    elif args.quiet:
        verbosidad = 0
    else:
        verbosidad = 1
    print('verbosidad', verbosidad)

    # Comprobamos si queremos/hace falta calcular los bias/flats
    print('bias:', args.nobias, args.sibias)
    realizarbias = decidir_repetir_calculos(args.nobias, args.sibias, 'bias', args.dir_dataf, args.dir_bias)
    print('flat:', args.noflat, args.siflat)
    realizarflat = decidir_repetir_calculos(args.noflat, args.siflat, 'flat', args.dir_dataf, args.dir_flats)

    # Creamos una lista de las noches disponibles
    lista_noches = os.listdir(args.dir_datos)
    lista_noches.sort()
    tiempo_inicio = time.time()

    # Separamos entre calibración y ciencia
    create_list_cal_and_sci(lista_noches, args.dir_listas, args.dir_datos, desc_bias, desc_flats, desc_arc,
                            verbosidad, args.calysci)
    tiempo_listas = time.time()

    print(args.nobias, args.noflat, args.noreducc)

    # importado_b = pd.read_csv('df_bias.csv')
    # importado_f = pd.read_csv('df_flat.csv')
    # pruebas_pandas(importado_f)

    # Creamos los Master Biases
    if realizarbias:
        df_bias = make_master_bias(lista_noches, args.dir_listas, args.dir_datos, args.dir_bias,
                                   args.interactive, args.recortar, verbose_imagen=args.verboseimage)
        numero_bias = len(os.listdir(args.dir_bias))
        print(df_bias)
        _ = df_bias.to_csv('df_bias.csv', index=None, header=True)
    else:
        df_bias = pd.read_csv('df_bias.csv')
        print('Se han importado los bias')
        numero_bias = '-'

    lista_bias = obtain_files_lists(args.dir_bias)
    tiempo_biases = time.time()

    # Creamos los Master Flats
    if realizarflat:
        df_flat = make_master_flat(lista_noches, lista_bias,
                                   args.dir_listas, args.dir_datos, args.dir_bias, args.dir_flats,
                                   verbosidad, args.interactive, verbose_imagen=args.verboseimage)
        numero_flats = len(os.listdir(args.dir_flats))
        print(df_flat)
        _ = df_flat.to_csv('df_flat.csv', index=None, header=True)
    else:
        df_flat = pd.read_csv('df_flat.csv')
        print('Se han importado los flats')
        numero_flats = '-'

    tiempo_flats = time.time()

    # Juntamos todos los procesos y relizamos la reducción
    if args.noreducc:
        numeros_reducidos = reducing_images(lista_noches, args.dir_listas, args.dir_datos, args.dir_bias,
                                            args.dir_flats, args.dir_reducc, df_bias, df_flat,
                                            verbosidad, verbose_imagen=args.verboseimage)
    else:
        numeros_reducidos = '-'
    tiempo_reducc = time.time()

    # Mostramos resultados de ambos procesos
    mostrarresultados(['Tiempo Listas', 'Tiempo Master Bias', 'Tiempo Master Flats', 'Tiempo Reduccion', 'Tiempo Total',
                       'Cuantos Biases', 'Cuantos Flats', 'Cuantos Reducidos'],
                      [round(tiempo_listas - tiempo_inicio, 2), round(tiempo_biases - tiempo_listas, 2),
                       round(tiempo_flats - tiempo_biases, 2), round(tiempo_reducc - tiempo_flats, 2),
                       round(tiempo_reducc - tiempo_listas, 2),
                       numero_bias, numero_flats, numeros_reducidos],
                      titulo='Tiempo Ejecucion')


if __name__ == "__main__":

    main()
