import json, os
m='e:/DWCT/OLD/Others/PD_Project_2026/Proteomics/PXD001676/dump/PXD001676_manifest.json'
with open(m) as f:
    records=json.load(f)
missing=[]
present=[]
for r in records:
    if r.get('category')=='SEARCH' and r['fileName'].endswith('.xml'):
        p=os.path.join('e:/DWCT/OLD/Others/PD_Project_2026/Proteomics/PXD001676', r['fileName'])
        if os.path.exists(p): present.append(r['fileName'])
        else:
            p2=os.path.join('e:/DWCT/OLD/Others/PD_Project_2026/Proteomics/PXD001676/dump', r['fileName'])
            if os.path.exists(p2): present.append(r['fileName'])
            else: missing.append(r['fileName'])
print('present:', len(present))
for x in sorted(present): print(x)
print('\nmissing:', len(missing))
for x in sorted(missing): print(x)
