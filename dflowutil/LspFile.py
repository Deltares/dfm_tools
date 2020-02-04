import pandas as pd

class LspFile():
    '''
    a water quality model *.lsp file
    '''

    def __init__(self, lspfile, procfile):
        """

        Arguments:
            lspfile {str} -- path to lsp file
            procfile {str} -- path to proces(m).asc file
        """
        self.path = lspfile
        self.proc = procfile
        self.get_units()

    def lsp_extract(self, var):
        """
        parses a line in an lsp file

        Arguments:
            var {str} -- line in an lsp file
        """
        ind = [ii for ii, jj in enumerate(var) if jj == '[' or jj == ']']
        name = var[ind[0] + 1:ind[1]].strip()
        descript = var[ind[1] + 1:-1].strip()
        return name, descript

    def get_units(self):

        with open(self.proc, 'r') as proc:
            unit = {}
            page = proc.readlines()
            for ind, line in enumerate(page):
                if line[0:10].strip() not in unit.keys():
                    print(line)
                    unit[line[0:10].strip()] = line[85:].strip()
                    print(line[85:].strip())
        self.units = unit

    def lsp_to_table(self, tablefile):
        """

        Arguments:
            tablefile {str} -- path to output csv file
        """
        self.tablefile = tablefile
        with open(self.path, 'r') as lsp:
            with open(self.tablefile, 'w') as table:
                table.write('process,process description,parameter,value,unit,description,\n')
                page = lsp.readlines()
                procInd = []
                sectName = []
                sectInd = []
                for ind, line in enumerate(page):
                    if '#' in line:
                        sectName.append(line)
                        sectInd.append(ind)
                sectInd0 = [ii for ii, jj in enumerate(page) if 'determining the input for the processes' in jj]
                for ind, line in enumerate(page):
                    if "Input for [" in line and ind > sectInd0[0]:
                        procInd.append(ind)
                for proc in range(0, len(procInd) - 1):
                    locProc = page[procInd[proc]:procInd[proc + 1]]
                    procName, procDescript = self.lsp_extract(locProc[0])
                    for ll in range(1, len(locProc) - 2, 2):
                        if ll == 1:
                            table.write(('%s,%s,') % (procName.replace(',', ''), procDescript.replace(',', '')))
                        else:
                            table.write(' , ,')
                        parName, parDescript = self.lsp_extract(locProc[ll])
                        table.write(('%s,') % (parName))
                        valLine = locProc[ll + 1]
                        if 'using output from' in valLine:
                            table.write(('%s') % valLine.replace(',', '').replace('\n', '').strip())
                        elif 'using substance' in valLine:
                            table.write(valLine.replace('Using ', '').replace(',', '').replace('\n', '').strip())
                        elif 'using constant nr' in valLine:
                            col = [ii for ii, jj in enumerate(valLine) if jj == ':']
                            table.write(('%s') % valLine[col[0] + 2:].replace(',', '').replace('\n', '').strip())
                        elif ':' in valLine:
                            col = [ii for ii, jj in enumerate(valLine) if jj == ':']
                            table.write(('%s') % valLine[col[0] + 2:].replace(',', '').replace('\n', '').strip())
                        else:
                            table.write(valLine.replace('\n', '').strip())
                        if parName != 'fcPPGreeN':
                            try:
                                table.write(',%s,%s\n' % (self.units[parName], parDescript.replace(',', '')))
                            except:
                                print('No unit found for %s in given process library' % parName)
                                table.write(',NOUNIT,%s\n' % (parDescript.replace(',', '')))
                        else:
                            table.write(',%s,%s\n' % ('(gC/m3/d)', parDescript.replace(',', '')))

    # create latex table
    def lsp_to_latex(self, latexfile):
        """

        Arguments:
            latexfile {str} -- path to output latex table file
        """
        self.latexfile = latexfile
        dat = pd.read_csv(self.tablefile)
        if isinstance(self.latexfile, str):
            dat.drop('Unnamed: 6', axis=1, inplace=True)
            bc = ['process', 'process description']
            with open(latexfile, 'w') as ltx:
                ltx.write('\\\begin{longtable}{')
                for cc in dat.columns:
                    if cc not in bc:
                        ltx.write('|l')
                ltx.write('|')
                ltx.write('}\n')
                for ii, cc in enumerate(dat.columns):
                    if cc not in bc and ii != len(dat.columns):
                        ltx.write('\\\textbf{' + cc + '} & ')
                    elif cc not in bc:
                        ltx.write('\\\textbf{' + cc + '} ')

                ltx.write("\\\\ \n")
                subs = {}
                for ii, ser in enumerate(dat[dat.columns[3]]):
                    sub = dat['parameter'].iloc[ii]
                    if sub not in subs.keys() and 'Using' not in dat['value'].iloc[ii]:
                        for ind, val in enumerate(dat.iloc[ii]):
                            if dat.columns[ind] not in bc:
                                if ind != len(dat.iloc[ii]) - 1:
                                    ltx.write(str(val).replace('_', '') + ' & ')
                                else:
                                    ltx.write(str(val).replace('_', ''))

                        subs[sub] = sub
                        ltx.write("\\\\ \n")
                ltx.write('\\end{longtable}')
