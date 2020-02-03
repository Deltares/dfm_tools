class OutFile():
    '''
    dflow outfile
    '''
    def __init__(self, path):
        """
        Arguments:
            path {str} -- path to out.txt file
        """
        self.path = path
        with open(path, 'r') as out:
            lines = out.readlines()
            self.warnings = [line for line in lines if 'WARNING' in line]
        self.sours_warning = [line for line in self.warnings if 'outside model area' in line]

