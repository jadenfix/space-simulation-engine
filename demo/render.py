"""Render actual periodic plasma benchmark convergence results, not particle trajectories.
Run from repo root: python3 demo/render.py (Pillow, make, gcc, ffmpeg, DejaVu fonts).
"""
import pathlib,hashlib,json,subprocess,tempfile,math
from PIL import Image,ImageDraw,ImageFont
ROOT=pathlib.Path(__file__).resolve().parents[1]
manifest=json.loads((ROOT/'demo/source-files.json').read_text())
for item in manifest['files']:
    b=(ROOT/item['path']).read_bytes()
    actual=hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()
    if actual!=item['sha']:raise SystemExit('Source changed: reverify provenance before rendering '+item['path'])
p=subprocess.run(['make','-C','assurance','-j2','test'],cwd=ROOT,capture_output=True,text=True,check=True)
(ROOT/'demo/core-run.txt').write_text(p.stdout+p.stderr)
receipt=json.loads((ROOT/'assurance/output/em_pic3d_plasma_modes_receipt.json').read_text())
assert receipt['passes'] and receipt['failures']==0 and all(c['passes'] for c in receipt['cases'].values())
(ROOT/'demo/plasma-results.json').write_text(json.dumps(receipt,indent=2)+'\n')
W,H,FPS,SECONDS=960,540,12,16
font='/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'
mono='/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf'
def f(n,m=False):return ImageFont.truetype(mono if m else font,n)
cases=receipt['cases'];names=['coarse','medium','fine'];colors=['#ffb878','#80b9ff','#79dcca']
errors=[abs(cases[n]['measured_omega_rad_s']/cases[n]['continuum_omega_rad_s']-1) for n in names]
with tempfile.TemporaryDirectory() as td:
    td=pathlib.Path(td)
    for i in range(FPS*SECONDS):
        t=i/FPS;stage=min(int(t//4),3)
        im=Image.new('RGB',(W,H),'#0b1020');d=ImageDraw.Draw(im)
        d.rounded_rectangle((30,24,930,512),radius=22,fill='#141d32',outline='#283953')
        d.text((58,44),'Space Simulation Engine / numerical assurance',font=f(18),fill='#8eacc8')
        titles=['A periodic plasma reference benchmark','Refine the timestep. Watch error shrink.','Conservation, checked at every step','Evidence from the actual C solver']
        d.text((58,87),titles[stage],font=f(29),fill='#f2f6ff')
        d.text((58,133),'Synthetic initial conditions · periodic 8 × 1 × 1 grid · no device claim',font=f(16),fill='#8eacc8')
        if stage==0:
            rows=[('6 benchmark cases','electron / proton, density, mass, pair plasma'),('C solver + strict gates','continuity, Gauss, energy, momentum'),('Fresh source execution','make -C assurance -j2 test')]
            for k,(title,desc) in enumerate(rows):
                y=195+k*79;d.rounded_rectangle((58,y,900,y+65),radius=10,fill='#202d46')
                d.text((79,y+10),title,font=f(22),fill='#79dcca');d.text((350,y+17),desc,font=f(15),fill='#c2cfe0')
        elif stage==1:
            d.text((60,180),'RELATIVE FREQUENCY ERROR VS CONTINUUM REFERENCE',font=f(15,True),fill='#8eacc8')
            reveal=min(3,1+int((t-4)//1.2))
            for k,n in enumerate(names[:reveal]):
                y=229+k*68;length=int(540*errors[k]/max(errors))
                d.text((60,y),n,font=f(18),fill=colors[k]);d.rounded_rectangle((165,y,165+length,y+25),radius=5,fill=colors[k])
                d.text((165+length+12,y),f'{errors[k]:.3e}',font=f(17,True),fill='#e4eefb')
            d.text((60,434),'Observed refinement order: '+f"{receipt['derived']['order_coarse_medium']:.4f} / {receipt['derived']['order_medium_fine']:.4f}",font=f(22),fill='#79dcca')
        elif stage==2:
            d.text((60,180),'PAIR-PLASMA CASE · MAXIMUM RELATIVE RESIDUALS',font=f(15,True),fill='#8eacc8')
            p=cases['pair_plasma']
            metrics=[('Total energy drift','maximum_total_energy_relative_drift'),('Charge continuity','maximum_continuity_relative'),('Gauss constraint','maximum_gauss_relative'),('Momentum closure','maximum_momentum_relative')]
            for k,(label,key) in enumerate(metrics):
                y=223+k*48;d.text((65,y),label,font=f(21),fill='#c2cfe0');d.text((575,y),f'{p[key]:.3e}',font=f(23,True),fill='#79dcca')
            d.text((65,435),'All six emitted benchmark cases: PASS',font=f(21),fill='#79dcca')
        else:
            d.text((60,180),'MEASURED FREQUENCY · RAD/S',font=f(16,True),fill='#8eacc8')
            for k,n in enumerate(names):
                y=225+k*53;d.text((65,y),n,font=f(22),fill=colors[k]);d.text((325,y),f"{cases[n]['measured_omega_rad_s']:.6e}",font=f(23,True),fill='#e4eefb')
            d.text((65,413),'Reference benchmark validated; propulsion remains unproven.',font=f(19),fill='#ffb878')
            d.text((65,447),'Generator + raw JSON receipt included in this repository.',font=f(17),fill='#c2cfe0')
        d.text((58,481),'Rendered actual convergence outputs · source '+manifest['source_commit'][:12],font=f(13),fill='#8eacc8')
        d.rounded_rectangle((30,529,30+int(900*(i+1)/(FPS*SECONDS)),533),radius=2,fill='#79dcca');im.save(td/f'{i:04}.png')
    subprocess.run(['ffmpeg','-v','error','-threads','1','-filter_threads','1','-filter_complex_threads','1','-y','-framerate',str(FPS),'-i',str(td/'%04d.png'),'-c:v','libx264','-crf','20','-pix_fmt','yuv420p','-movflags','+faststart',str(ROOT/'demo/plasma-convergence.mp4')],check=True)
    subprocess.run(['ffmpeg','-v','error','-threads','1','-filter_threads','1','-y','-framerate',str(FPS),'-i',str(td/'%04d.png'),'-vf','palettegen','-frames:v','1',str(td/'palette.png')],check=True)
    subprocess.run(['ffmpeg','-v','error','-threads','1','-filter_complex_threads','1','-y','-framerate',str(FPS),'-i',str(td/'%04d.png'),'-i',str(td/'palette.png'),'-lavfi','paletteuse','-loop','0',str(ROOT/'demo/plasma-convergence.gif')],check=True)
    sheet=Image.new('RGB',(960,540))
    for k,n in enumerate([24,80,120,170]):sheet.paste(Image.open(td/f'{n:04}.png').resize((480,270)),((k%2)*480,(k//2)*270))
    sheet.save(ROOT/'demo/contact.png')
print('PASS: source manifest, real suite, all six benchmark receipts; GIF and MP4 generated')
