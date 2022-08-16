import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("csv", help="Input csv", type=str)
    return parser.parse_args()


def main():
    args = parse_args()
    output = ""
    header_flag = False
    with open(args.csv, "r") as f1:
        for line in f1:
            line_split = line.split(",")
            if not header_flag:
                output += "\\begin{table*}[t]\n\\centering\n\\scriptsize\n\\begin{tabular}{"
                output += "|c" * len(line_split)
                output += "|}\n\\hline\n"
                for i in line_split:
                    output += "\\textbf{" + i + "} &"
                output = output[:-2]
                output += "\\\\\n\\hline\n\\hline\n\n"
                header_flag = True
            else:
                for i in line_split:
                    output += i + " & "
                output = output[:-2]
                output += "\\\\\n\\hline\n"
        output += "\n\n\\hline\n\\end{tabular}\n\\label{tab1}\n\\caption{\\textbf{TEXT HERE} TEXT HERE}\n"
        output += "\\end{table*}"
    with open("output_latex.txt", "w") as f2:
        f2.write(output)


if __name__ == "__main__":
    main()
