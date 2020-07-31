"""Summary
"""


def read_file(name_file):
    """Take in list file and output decoy/conformer name list.

    Args:
        name_file (str): name of file

    Returns:
        list: list of row in the input file
    """
    name_list = []
    with open(name_file, "r") as open_file:
        for line in open_file:
            name = line.strip()
            name_list.append(name)

        open_file.close()

    return name_list
