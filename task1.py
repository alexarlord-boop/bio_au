from Utils import IOUtils
import string


def format_strings(data):
    dm = len(data)
    for i in range(dm):
        ds = len(data[i])
        line = list(data[i])
        for j in range(ds):
            if line[j] in string.ascii_letters or line[j] in string.digits or line[j] in [">", "_", "\n"]:
                continue
            else:
                line[j] = "_"

        data[i] = "".join(line)
    return data


def cut_strings(data):
    dm = len(data)
    for i in range(dm):
        if data[i][0] == ">":
            line = list(data[i])
            if len(line) > 99:
                data[i] = line[0:98:] + ["\n"]
    return data


def combine_target_and_query(target, query):
    return target[:-1] + query


if __name__ == '__main__':
    query_in = 'C:\\bio\\query_in.txt'
    query_out = 'C:\\bio\\query_out.txt'
    target_in = "C:\\bio\\target_in.txt"
    target_out = "C:\\bio\\target_out.txt"
    total_out = "C:\\bio\\total.txt"
    io = IOUtils()

    qdata = io.open_file(query_in)
    tdata = io.open_file(target_in)
    io.write_data(target_out, format_strings(cut_strings(tdata)))
    io.write_data(query_out, format_strings(cut_strings(qdata)))
    io.write_data(total_out, combine_target_and_query(tdata, qdata))
