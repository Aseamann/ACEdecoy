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
            line_split = line[:-1].split(",")
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
                color_flag = False
                for i in range(len(line_split)):
                    if not color_flag:
                        output += "\\textcolor{"
                        if float(line_split[i + 1]) <= -1.0:
                            output += "red}{" + line_split[i]
                        else:
                            output += "orange}{" + line_split[i]
                        output += "} & "
                        color_flag = True
                    else:
                        output += line_split[i] + " & "
                output = output[:-3]
                output += "\\\\\n\\hline\n"
        output += "\\end{tabular}\n\\label{tab1}\n\\caption{\\textbf{TEXT HERE} TEXT HERE}\n"
        output += "\\end{table*}"
    with open("output_latex.txt", "w") as f2:
        f2.write(output)


if __name__ == "__main__":
    main()
