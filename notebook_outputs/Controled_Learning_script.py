import csv
import numpy as np
import ast
import re
import matplotlib.pyplot as plt
import pandas as pd

with open('tissue_cell_pairs.tsv', newline='') as file:
    reader = csv.reader(file, delimiter='\t')
    tissue_cells = list(reader)
tissue_cells = tissue_cells[1:]
print("number of cleaned cells:", len(tissue_cells))

with open('common_cells_across_tissues.csv', newline='') as file:
    reader = csv.reader(file, delimiter='\t')
    common_cells = list(reader)
common_cells.append('serous glandular cells')
print("number of common cells:", len(common_cells))

with open('gene_list.csv', newline='') as file:
    reader = csv.reader(file, delimiter='\t')
    gene_list = list(reader)
print("number of cleaned genes:", len(gene_list))

# --- Load nTPM matrices separately ---

# Load high method nTPM matrix
nTPM_matrix_high_df = pd.read_csv('gene_expression_matrix_high.csv', sep='\t')
nTPM_matrix_high_only = nTPM_matrix_high_df.drop(columns=["Tissue", "Cell type"])
nTPM_matrix_high = nTPM_matrix_high_only.values

# Load median method nTPM matrix
nTPM_matrix_median_df = pd.read_csv('gene_expression_matrix_median.csv', sep='\t')
nTPM_matrix_median_only = nTPM_matrix_median_df.drop(columns=["Tissue", "Cell type"])
nTPM_matrix_median = nTPM_matrix_median_only.values

# --- Precompute mappings ---
if isinstance(gene_list[0], list):
    gene_list = [gene[0] for gene in gene_list]
tissue_cell_to_index = {tuple(row): idx for idx, row in enumerate(tissue_cells)}
gene_to_index = {gene: idx for idx, gene in enumerate(gene_list)}

# --- Load positive labels ---
positives_df = pd.read_csv('positives_labels.csv', sep='\t')
print(f"Number of cells with positive labels: {len(positives_df)}")
positive_label_count = positives_df["Positive Gene Names"].apply(lambda x: len(eval(x))).sum()
print(f"Total number of positive labels: {positive_label_count}")
# Formulate positive labels using indexes
positives= {}
for idx, row in positives_df.iterrows():
    tissue = row["Tissue"]
    cell_type = row["Cell type"]
    gene_names = eval(row["Positive Gene Names"])

    key = (tissue, cell_type)
    if key not in tissue_cell_to_index:
        continue

    cell_idx = tissue_cell_to_index[key]
    gene_indices = [gene_to_index[gene] for gene in gene_names if gene in gene_to_index]

    if gene_indices:
        positives[cell_idx] = gene_indices
print(f"Number of cells in positives (indexed): {len(positives)}")


# --- Load negative labels ---
negatives_df = pd.read_csv('negative_labels.csv', sep='\t')
print(f"Number of cells with negative labels: {len(negatives_df)}")
negative_label_count = negatives_df["Negative Gene Names"].apply(lambda x: len(eval(x))).sum()
print(f"Total number of negative labels from list: {negative_label_count}")
# Formulate negative labels using indexes
negatives = {}
for idx, row in negatives_df.iterrows():
    tissue = row["Tissue"]
    cell_type = row["Cell type"]
    gene_names = eval(row["Negative Gene Names"])

    key = (tissue, cell_type)
    if key not in tissue_cell_to_index:
        continue

    cell_idx = tissue_cell_to_index[key]
    gene_indices = [gene_to_index[gene] for gene in gene_names if gene in gene_to_index]

    if gene_indices:
        negatives[cell_idx] = gene_indices
print(f"Number of cells in negatives (indexed): {len(negatives)}")


nTPM_matrix_high_medians = np.median(nTPM_matrix_high, axis=0).tolist()
genes_sorted_high_ = [index for index, value in sorted(enumerate(nTPM_matrix_high_medians ), key=lambda x: x[1], reverse=True)]
negatives1_high = genes_sorted_high_[:50]

nTPM_matrix_median_medians = np.median(nTPM_matrix_high, axis=0).tolist()
genes_sorted_median = [index for index, value in sorted(enumerate(nTPM_matrix_high_medians ), key=lambda x: x[1], reverse=True)]
negatives1_median = genes_sorted_median[:50]

positives2 = [
    ["PECAM1", "lung", 'endothelial cells', "LNP", "PECAM-1 directed re-targeting of exogenous mRNA providing two orders of magnitude enhancement of vascular delivery and expression in lungs independent of apolipoprotein E-mediated uptake"],
    ["VCAM1", "vascular", 'endothelial cells', "LNP", "Selective targeting of nanomedicine to inflamed cerebral vasculature to enhance the blood–brain barrier"],
    ["CD4", "pbmc", "t-cells", "LNP", "Highly efficient CD4+ T cell targeting and genetic recombination using engineered CD4+ cell-homing mRNA-LNPs"],
    ["CD5", "pbmc", "t-cells", "LNP", "CAR T cells produced in vivo to treat cardiac injury"],
    ["CD19", "pbmc", "b-cells", "LNP", ""],
    ["CD3", "pbmc", "t-cells", "LNP", "conference; doudna paper"],
    ["NCR1", "pbmc", "nk-cells", "LNP", "conference"],
    ["CD14", "pbmc", "macrophages", "LNP", "conference"],
    ["MRC1", "pbmc", "macrophages", "LNP", "conference"],
    ["ITGAM", "pbmc", "macrophages", "", "bacteria injector"],
    ["CD28", "pbmc", "t-cells", "EDV", "doudna EDV"],
    ["CD40", "pbmc", "b-cells", "lenti", "fengzhang new lenti papr"],
    ["ENG", "vascular", "endothelial cells", "LNP", "Targeting of immunoliposomes to endothelial cells using a single-chain Fv fragment directed against human endoglin (CD105)"],
    ["MRC1", "pbmc", "dendritic cells", "LNP", "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9322927/#B8-pharmaceuticals-15-00897"],
    ["CD8", "pbmc", "t-cells", "LNP", ""],
    ["PDPN", "skin", "endothelial cells", "LNP", "Targeted delivery of lipid nanoparticle to lymphatic endothelial cells via anti-podoplanin antibody"],
    ["PLVAP", "lung", "endothelial cells", "LNP?", "https://pubs.acs.org/doi/full/10.1021/acschembio.0c00003"],
    ["FCER2", "pbmc", "b-cells", "", "WildDISCO whole body imaging"]
]

# nTPM_matrix_high
# Train for parameters p, q, r

m = len(tissue_cells)
dic = {}
min_score = 0
for p in np.linspace(-700, -500, 20):
  for q in np.linspace(-700, -500, 20):
    for r in np.linspace(-700, -20, 20):
      penalty_matrix = [] #penalty for jth tissue_cell based on ith tissue_cell
      for i in range(m):
        row = []
        for j in range(m):
          if j == i:
            row.append(1000)
          elif tissue_cells[j][1] == tissue_cells[i][1]:
            if tissue_cells[j][1] in common_cells:
              row.append(0) # same cell (common) and diff tissue
            else:
              row.append(p) # same cell (non common) and diff tissue
          elif tissue_cells[j][0] == tissue_cells[i][0]:
            row.append(q) # diff cell and same tissue
          else:
            row.append(r) # diff cell and diff tissue
        penalty_matrix.append(row)
      penalty_matrix = np.array(penalty_matrix)
      objective_matrix = np.dot(penalty_matrix, nTPM_matrix_high)
      objective_matrix = objective_matrix.tolist()

      count1 = 0
      for cell in positives:
        markers = positives[cell]
        obj_row = objective_matrix[cell]
        markers_suggested = sorted(list(enumerate(obj_row)), key = lambda x: x[1], reverse=True)[:10]
        markers_suggested = [index for index, value in markers_suggested]
        count1 += len(set(markers).intersection(set(markers_suggested)))

      count2 = 0
      for row in positives2:
        if [row[0]] in gene_list:
          marker = gene_list.index([row[0]])
          cell = tissue_cells.index([row[1], row[2]])
          obj_row = objective_matrix[cell]
          markers_suggested = sorted(list(enumerate(obj_row)), key = lambda x: x[1], reverse=True)[:10]
          markers_suggested = [index for index, value in markers_suggested]
          if marker in markers:
            count2 += 1

      count3 = 0
      for cell in negatives:
        non_markers = negatives[cell]
        obj_row = objective_matrix[cell]
        markers_suggested = sorted(list(enumerate(obj_row)), key = lambda x: x[1], reverse=True)[:10]
        markers_suggested = [index for index, value in markers_suggested]
        count3 += len(set(non_markers).intersection(set(markers_suggested)))

      count4 = 0
      for obj_row in objective_matrix:
        markers_suggested = sorted(list(enumerate(obj_row)), key = lambda x: x[1], reverse=True)[: 10]
        markers_suggested = [index for index, value in markers_suggested]
      for i in markers_suggested:
        if i in negatives1_high:
          count4 += 1


      if count1 + count2 - count3 - count4 > min_score:
        min_score = count1 + count2 - count3 - count4
        max_set = (p, q, r)

      dic[(p, q, r)] = (count1, count2, count3, count4)

print("The (nearly-) optimal parameters p, q, r:", max_set)
print("The number of hitted top 50 highly expression genes:", count4)
print("The number of hitted non- markers:", count3)
print("The number of hitted whole-body markers:", count2)
print("The number of hitted organ-wide markers:", count1)

# Form the score table
p = max_set[0]
q = max_set[1]
r = max_set[2]
m = len(tissue_cells)
t = 1
penalty_matrix = [] #penalty for jth tissue_cell based on ith tissue_cell
for i in range(m):
  row = []
  for j in range(m):
    if j == i:
      row.append(1000)
    elif tissue_cells[j][1] == tissue_cells[i][1]:
      if tissue_cells[j][1] in common_cells:
        row.append(0) # same cell (common) and diff tissue
      else:
        row.append(p) # same cell (non common) and diff tissue
    elif tissue_cells[j][0] == tissue_cells[i][0]:
      row.append(q) # diff cell and same tissue
    else:
      row.append(r) # diff cell and diff tissue
  penalty_matrix.append(row)
penalty_matrix = np.array(penalty_matrix)
objective_matrix = np.dot(penalty_matrix, nTPM_matrix_high)

# Return the top 10 recommeded markers for all cells
whole_body_markers = {}
for i in range(m):
  obj_row = objective_matrix[i]
  markers = sorted(list(enumerate(obj_row)), key = lambda x: x[1], reverse=True)[: 10]
  markers_id = [x[0] for x in markers]
  print()
  whole_body_markers[tissue_cells[i][0] + " " + tissue_cells[i][1]] = [genes[x][0] for x in markers_id]

filename = 'recommended_whole_body_markers_high.csv'
with open(filename, 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(["cell", "markers"])
    for key, value in whole_body_markers.items():
        row = [key, value]
        csvwriter.writerow(row)

# nTPM_matrix_median
# Train for parameters p, q, r
m = len(tissue_cells)
dic = {}
min_score = 0
for p in np.linspace(-700, -500, 20):
  for q in np.linspace(-700, -500, 20):
    for r in np.linspace(-700, -20, 20):
      penalty_matrix = [] #penalty for jth tissue_cell based on ith tissue_cell
      for i in range(m):
        row = []
        for j in range(m):
          if j == i:
            row.append(1000)
          elif tissue_cells[j][1] == tissue_cells[i][1]:
            if tissue_cells[j][1] in common_cells:
              row.append(0) # same cell (common) and diff tissue
            else:
              row.append(p) # same cell (non common) and diff tissue
          elif tissue_cells[j][0] == tissue_cells[i][0]:
            row.append(q) # diff cell and same tissue
          else:
            row.append(r) # diff cell and diff tissue
        penalty_matrix.append(row)
      penalty_matrix = np.array(penalty_matrix)
      objective_matrix = np.dot(penalty_matrix, nTPM_matrix_median)
      objective_matrix = objective_matrix.tolist()

      count1 = 0
      for cell in positives:
        markers = positives[cell]
        obj_row = objective_matrix[cell]
        markers_suggested = sorted(list(enumerate(obj_row)), key = lambda x: x[1], reverse=True)[:10]
        markers_suggested = [index for index, value in markers_suggested]
        count1 += len(set(markers).intersection(set(markers_suggested)))

      count2 = 0
      for row in positives2:
        if [row[0]] in gene_list:
          marker = gene_list.index([row[0]])
          cell = tissue_cells.index([row[1], row[2]])
          obj_row = objective_matrix[cell]
          markers_suggested = sorted(list(enumerate(obj_row)), key = lambda x: x[1], reverse=True)[:10]
          markers_suggested = [index for index, value in markers_suggested]
          if marker in markers:
            count2 += 1

      count3 = 0
      for cell in negatives:
        non_markers = negatives[cell]
        obj_row = objective_matrix[cell]
        markers_suggested = sorted(list(enumerate(obj_row)), key = lambda x: x[1], reverse=True)[:10]
        markers_suggested = [index for index, value in markers_suggested]
        count3 += len(set(non_markers).intersection(set(markers_suggested)))

      count4 = 0
      for obj_row in objective_matrix:
        markers_suggested = sorted(list(enumerate(obj_row)), key = lambda x: x[1], reverse=True)[: 10]
        markers_suggested = [index for index, value in markers_suggested]
      for i in markers_suggested:
        if i in negatives1_median:
          count4 += 1


      if count1 + count2 - count3 - count4 > min_score:
        min_score = count1 + count2 - count3 - count4
        max_set = (p, q, r)

      dic[(p, q, r)] = (count1, count2, count3, count4)

print("The (nearly-) optimal parameters p, q, r:", max_set)
print("The number of hitted top 50 highly expression genes:", count4)
print("The number of hitted non- markers:", count3)
print("The number of hitted whole-body markers:", count2)
print("The number of hitted organ-wide markers:", count1)

# Form the score table
p = max_set[0]
q = max_set[1]
r = max_set[2]
m = len(tissue_cells)
t = 1
penalty_matrix = [] #penalty for jth tissue_cell based on ith tissue_cell
for i in range(m):
  row = []
  for j in range(m):
    if j == i:
      row.append(1000)
    elif tissue_cells[j][1] == tissue_cells[i][1]:
      if tissue_cells[j][1] in common_cells:
        row.append(0) # same cell (common) and diff tissue
      else:
        row.append(p) # same cell (non common) and diff tissue
    elif tissue_cells[j][0] == tissue_cells[i][0]:
      row.append(q) # diff cell and same tissue
    else:
      row.append(r) # diff cell and diff tissue
  penalty_matrix.append(row)
penalty_matrix = np.array(penalty_matrix)
objective_matrix = np.dot(penalty_matrix, nTPM_matrix_median)

# Return the top 10 recommeded markers for all cells
whole_body_markers = {}
for i in range(m):
  obj_row = objective_matrix[i]
  markers = sorted(list(enumerate(obj_row)), key = lambda x: x[1], reverse=True)[: 10]
  markers_id = [x[0] for x in markers]
  print()
  whole_body_markers[tissue_cells[i][0] + " " + tissue_cells[i][1]] = [genes[x][0] for x in markers_id]

filename = 'recommended_whole_body_markers_median.csv'
with open(filename, 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(["cell", "markers"])
    for key, value in whole_body_markers.items():
        row = [key, value]
        csvwriter.writerow(row)

# Screen the cells and associted recommended markers meet the condition1: 1.5x
cells1 = set()
for i in range(m):
  top_3_indices = objective_matrix[i,:].argsort()[-3:][::-1]
  for j in top_3_indices:
    cell = np.argmax(nTPM_matrix[:, j])
    if cell == i:
      column_data = nTPM_matrix[:, j].copy()
      column_data[i] = -np.inf
      off_target = np.argmax(column_data)
      if nTPM_matrix[i, j] > 1.5 * nTPM_matrix[i, off_target]:
        cells1.add(i)

# Screen the cells and associted recommended markers meet the condition2: 1.5x
cells2 = set()
for i in range(m):
  top_3_indices = objective_matrix[i,:].argsort()[-3:][::-1]
  for j in top_3_indices:
    column_data = nTPM_matrix[:, j].copy()
    column_data = np.delete(column_data, i)
    average_value = np.mean(column_data)
    if nTPM_matrix[i, j] > 10 * average_value:
      cells2.add(i)

# Combine
union = cells1.union(cells2)

print("The number of cells satisfies condition1:", len(cells1))
print("The number of cells satisfies condition2:", len(cells2))
print("The number of cells in union set:", len(union))

# Return the top 10 recommeded markers for selected cells
whole_body_markers = {}
for i in union:
  obj_row = objective_matrix[i]
  markers = sorted(list(enumerate(obj_row)), key = lambda x: x[1], reverse=True)[: 10]
  markers_id = [x[0] for x in markers]
  print()
  whole_body_markers[tissue_cells[i][0] + " " + tissue_cells[i][1]] = [genes[x][0] for x in markers_id]

filename = 'selected_cells_recommended_whole_body_markers.csv'
with open(filename, 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(["cell", "markers"])
    for key, value in whole_body_markers.items():
        row = [key, value]
        csvwriter.writerow(row)



positives_df
