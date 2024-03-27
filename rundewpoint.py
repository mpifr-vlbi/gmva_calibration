#!/usr/bin/env python3
import re
import math
import subprocess
import os

def format_time_component(component):
    return component.zfill(2)

def calculate_dewpoint(temp, humi):
    ES = 6.11 * 10 ** (7.5 * temp / (237.7 + temp))
    E = humi * ES / 100.0
    y = math.log(E / 6.11) / 7.5 / math.log(10)
    dewtemp = 237.7 * y / (1 - y)
    return dewtemp

def dewpoint(infile, statcode):
    subprocess.run(['grep', '\/wx\/', infile], stdout=open('dum', 'w'))

    with open('dum', 'r') as f:
        formatted_content = [line.replace('wx/', 'wx/ ').replace(',', ' ') for line in f.readlines()]
    with open(infile + '.wx.form', 'w') as f:
        f.write(''.join(formatted_content))

    wx_data = []
    with open(infile + '.wx.form', 'r') as f_in:
        wx_lines = f_in.readlines()
        for line in wx_lines:
            parts = line.split()
            date_time = parts[0].split('.')
            day = format_time_component(date_time[1])
            hh = format_time_component(date_time[2][:2])
            mm = format_time_component(date_time[2][3:5])
            ss = format_time_component(date_time[2][6:8])
            temp = float(parts[1])
            press = float(parts[2])
            humi = float(parts[3])

            dewtemp = calculate_dewpoint(temp, humi)

            wx_data.append(f"{day}-{hh}:{mm}:{ss}   {temp:.1f}   {press:.1f}    {dewtemp:.1f}  0.0  0.0  0.00  0.0\n")

    with open('WX.' + statcode, 'w') as f:
        f.write(f"WEATHER {statcode} /\n")
        f.writelines(wx_data)
        f.write(" / \n")
    print(f"Weather file in WX.{statcode}")
    subprocess.run(['head', 'WX.' + statcode])
    os.remove('dum')
    os.remove(infile + '.wx.form')



def dewpoint_vlba(infile):
    subprocess.run(['awk', '{print > out}; /! The file contains weather data collected at the NRAO stations./{out="vlba_weather.txt"}', f'out={infile.replace(".vlba", "_tsys.txt")}', infile], check=True)

    stations = ["BR", "FD", "KP", "LA", "MK", "NL", "OV", "PT"]

    for station in stations:
        results = subprocess.run(['grep', '-n', station, 'vlba_weather.txt'], stdout=subprocess.PIPE, check=True, text=True)
        results_line = results.stdout.splitlines()[0].split(':')[0]

        with open('vlba_weather.txt', 'r') as f_in, open('temp.txt', 'w') as f_out:
            for line_num, line in enumerate(f_in, start=1):
                if line_num == int(results_line):
                    f_out.write(f"/ \nWEATHER {station} /\n")
                f_out.write(line)

        os.replace('temp.txt', 'vlba_weather.txt')

    subprocess.run(['awk', '/^[^!]/ { print $0 }', 'vlba_weather.txt'], stdout=open('temp.txt', 'w'), check=True)
    os.replace('temp.txt', 'vlba_weather.txt')

    subprocess.run(['ed', 'vlba_weather.txt'], input=b'1d\nwq\n', check=True)
    subprocess.run(['sed', '-e', '$a\\', 'vlba_weather.txt'], stdout=open('temp.txt', 'w'), check=True)
    os.replace('temp.txt', 'vlba_weather.txt')
    with open('vlba_weather.txt', 'a') as f:
        f.write(' /\n')
    os.remove(infile.replace(".vlba", "_tsys.txt"))


nonTsysStar = ['EF', 'PV', 'MH']

for filename in os.listdir():
    if filename.endswith('cal.vlba'):
        dewpoint_vlba(filename)
    elif filename.endswith('.log'):
        stat_pattern =r'^.{5}([a-z]{2})' 
        match_stat = re.search(stat_pattern, filename)
        print(match_stat)
        if match_stat:
            station_code = match_stat.group(1).upper()
            print(station_code)
            if station_code in nonTsysStar:
                dewpoint(filename, station_code)

with open('SUM.WX', 'w') as f_sum:
    with open('vlba_weather.txt', 'r') as f_vlba:
        f_sum.write(f_vlba.read())
    for wx_file in os.listdir():
        if wx_file.startswith('WX.'):
            with open(wx_file, 'r') as f_wx:
                f_sum.write(f_wx.read())
