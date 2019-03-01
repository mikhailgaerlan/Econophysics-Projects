#!/usr/bin/env python
import csv
import dateutil
import dateutil.parser
import urllib

#ETF = 'SPY'
#filename = ETF + '_File.txt'

def get_etf(filename, ETF):
   data = {}
   #------------------------------------------------------------------#
   # Start with data from nasdaq -------------------------------------#
   base_url = 'http://ichart.finance.yahoo.com/table.csv?'  +\
                    's=' + ETF + '&a=0&b=1&c=1990'
       #
       # Open the stream and munch the header
   data_stream = urllib.urlopen(base_url)
   data_stream.next()
       #
       # Parse the data
   for row in csv.reader(data_stream):
       date = dateutil.parser.parse(row[0])
       key = date.strftime('%Y/%m/%d')
#       print row[:]
       if key not in data:
           if '-' not in row[1:5] and 0 not in map(float,row[1:5]):
               row[0] = date.year+int(date.strftime('%j'),10)/365.0
               row[5] = float(row[5])
               data[key] = str(row[0]) + '\t' + str(row[1]) + '\t' +str(row[2]) + '\t' +str(row[3]) + \
                            '\t' +str(row[4]) + '\t' +str(row[5]) + '\t' + str(key) + '\t' +str(row[6])
    #------------------------------------------------------------------#
   # Write output to file --------------------------------------------#
   with open(filename, 'w') as output:
       for key in sorted(data.keys()):
           output.write("%s\n" % data[key])

   #------------------------------------------------------------------#
   return None

if __name__ == "__main__":
   get_etf(filename, ETF)
