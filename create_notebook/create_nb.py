import nbformat as nbf

def create_notebook_cells(code_cells,notebook_filename):
    nb = nbf.v4.new_notebook()
    for code in code_cells:
        nb.cells.append(nbf.v4.new_code_cell(code))
    with open(notebook_filename, 'w') as f:
        nbf.write(nb, f)

