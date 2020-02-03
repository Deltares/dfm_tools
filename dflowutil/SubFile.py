class SubFile():
    '''
    delwaq substance file class
    '''

    def __init__(self, path):
        """
        Arguments:
            path {str} -- path to *.sub file
        """
        self.path = path

        with open(self.path, 'r') as subs:
            sub = []
            transportable = {}
            lines = subs.readlines()
            for line in lines:
                if line[0:9] == 'substance':
                    tmp = line.split(' ')
                    name = tmp[1].replace("'", '')
                    sub.append(name)
                    status = tmp[2].replace("'", '').replace('\n', '')
                    transportable[name] = status

        self.substances = sub
        self.transportable = transportable
