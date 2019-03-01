#  This code downloads historical price data from Yahoo! Finance
#
#  To download data sets, look up the name for the stock, index,
#  or ETF and insert in the file Stocks_List.txt
#
#  To use this code, change the file extension from .txt to .py
#
#  Then at the command line, enter: 
#       python Download_Stock_Data.py
#
#  Note that the columns are:
#
#  Decimal Date; Open; High; Low; Close; Volume; Analog Date; Adjusted Close
#
#!/usr/bin/env python

import sys 
import getopt
import write_stock_data
                                        
def main(argv=None):

    data_file = open('Stocks_List.txt', 'r')

    j_file=0
    for line in data_file:
        j_file += 1

    data_file.close()

    data_file = open('Stocks_List.txt', 'r')

    print ''
    print 'Stock or ETF Data Downloads: '

    j=0
    for line in data_file:
        j+=1
        items = line.strip().split(',')
    #
    #   percent_complete = 100.0*float(j)/float(j_file) 
    #   print '     Percent File Completed: ', "{0:.2f}".format(percent_complete),'%', "                                 \r",   

    #
    #   Both variables below are string variables
    #

        etf_name                = items[0]

        etf_output_file = etf_name + '_File.txt'

        print etf_name, '(#', j, 'of ', j_file, ')'

#
        write_stock_data.get_etf(etf_output_file,etf_name)
#
    data_file.close()
    print ''
#

if __name__ == "__main__": 
	sys.exit(main())


