class IOUtils:

    def __init__(self):
        pass

    def open_file(self, filename):
        data = []
        with open(filename) as f:
            data = f.readlines()
        return data

    def write_data(self, filename, data):
        with open(filename, 'w') as f:
            f.writelines(data)
            print("writing data to:" + filename)
