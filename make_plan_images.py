from glob import glob
from PIL import Image, ImageFont, ImageDraw
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

font_size = 12
font = 'data/FiraCode-Regular.ttf'
font = ImageFont.truetype(font, font_size)
w, h_ = 300, 300
vpadding = 8
h = h_ + (vpadding*2)+font_size

arr_w, arr_h = 48, 72
arr_thick = 4
arr_pad = 24
arr_head_h = 8
arr_head_w = 22
arr_center = arr_w/2-arr_thick/2
arrow = Image.new('RGB', (arr_w, arr_h), color='white')
draw = ImageDraw.Draw(arrow)
draw.line([(arr_center, 0+arr_pad), (arr_center, arr_h-arr_pad-arr_head_h)], width=arr_thick, fill='black')
draw.polygon([
    (arr_center-arr_head_w/2, arr_h-arr_pad-arr_head_h),
    (arr_center+arr_head_w/2, arr_h-arr_pad-arr_head_h),
    (arr_center, arr_h-arr_pad)], fill='black')


def make_image(smi, base=False):
    mol = Chem.MolFromSmiles(smi)
    formula = CalcMolFormula(mol)
    if base: formula = '{} (base)'.format(formula)
    mol_im = Draw.MolToImage(mol, size=(w, h_))
    im = Image.new('RGB', (w, h), color='white')
    im.paste(mol_im)
    draw = ImageDraw.Draw(im)
    tw, th = draw.textsize(formula)
    draw.text((w/2-tw/2, h-th-vpadding), formula, font=font, fill='black')
    return im

if __name__ == '__main__':
    import os
    import json
    from tqdm import tqdm
    from multiprocessing import Pool

    # Get most recent batch
    out_dir = sorted(os.listdir('data/sample'))[-1]
    files = glob('data/sample/{}/*.json'.format(out_dir))
    plans_dir = os.path.join('data/sample', out_dir, 'plans')
    if not os.path.exists(plans_dir):
        os.makedirs(plans_dir)

    def process(f):
        mols = json.load(open(f))
        for mol in mols:
            id = mol['image'].split('/')[-1].replace('.png', '')
            out_path = os.path.join(plans_dir, '{}.png'.format(id))
            if os.path.exists(out_path): continue
            if 'synthesis' not in mol: continue
            compounds = {}
            stages = []
            max_w = 0
            for step in mol['synthesis']['plan']:
                ims = []
                for m in step:
                    smi = m['smiles']
                    if smi not in compounds:
                        compounds[smi] = make_image(smi, m['is_terminal'])
                    ims.append(compounds[smi])

                tot_w = len(ims)*w
                img = Image.new('RGB', (tot_w, h), color='white')
                for i, im in enumerate(ims):
                    img.paste(im.copy(), (i*w, 0))
                if tot_w > max_w: max_w = tot_w
                stages.append(img)

            tot_h = h*len(stages) + arr_h*(len(stages)-1) + arr_h # last arr_h for extra bottom padding
            img = Image.new('RGB', (max_w, tot_h), color='white')
            cur_h = 0
            for i, im in enumerate(stages):
                if i > 0:
                    img.paste(arrow.copy(), (int(max_w/2 - arr_w/2), cur_h))
                    cur_h += arr_h
                img.paste(im, (int(max_w/2 - im.width/2), cur_h))
                cur_h += h

            img.save(out_path)

    with Pool() as p:
        for _ in tqdm(p.imap(process, files), total=len(files)): pass