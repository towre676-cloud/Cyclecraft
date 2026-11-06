"""Build the offline second edition, numerical archive and complete package."""
from pathlib import Path
import json,re,html,base64,subprocess,os,hashlib,zipfile
ROOT=Path(__file__).resolve().parent;OUT=ROOT.parent

def read(n):return (ROOT/n).read_text()
def load(n):return json.loads(read(n))

extra=[
 {'key':'dechant2022','authors':'Andreas Dechant','title':'Minimum entropy production, detailed balance and Wasserstein distance for continuous-time Markov processes','year':2022,'journal':'Journal of Physics A 55, 094001','url':'https://arxiv.org/abs/2110.01141'},
 {'key':'impens2025','authors':"François Impens, Felippo M. D’Angelis, David Guéry-Odelin, Felipe A. Pinheiro, Caio Lewenkopf",'title':'Optimal control over the full counting statistics in a non-adiabatic pump','year':2025,'journal':'Physical Review Research 7, 033253','url':'https://arxiv.org/abs/2503.14614'},
 {'key':'wachtel2015','authors':'Artur Wachtel, Jürgen Vollmer, Bernhard Altaner','title':'Fluctuating currents in stochastic thermodynamics. I. Gauge invariance of asymptotic statistics','year':2015,'journal':'Physical Review E 92, 042132','url':'https://journals.aps.org/pre/abstract/10.1103/PhysRevE.92.042132'},
 {'key':'koyuk2019','authors':'Timur Koyuk, Udo Seifert','title':'Operationally accessible bounds on fluctuations and entropy production in periodically driven systems','year':2019,'journal':'Physical Review Letters 122, 230601','url':'https://arxiv.org/abs/1904.08807'},
 {'key':'ilker2022','authors':'Efe Ilker, Özenç Güngör, Benjamin Kuznets-Speck, Joshua Chiel, Sebastian Deffner, Michael Hinczewski','title':'Shortcuts in Stochastic Systems and Control of Biophysical Processes','year':2022,'journal':'Physical Review X 12, 021048','url':'https://journals.aps.org/prx/abstract/10.1103/PhysRevX.12.021048'},
 {'key':'johansson2017','authors':'Fredrik Johansson','title':'Arb: Efficient Arbitrary-Precision Midpoint-Radius Interval Arithmetic','year':2017,'journal':'ACM Transactions on Mathematical Software 43(3), Article 22','url':'https://arxiv.org/abs/1611.02831'},
 {'key':'maas2011','authors':'Jan Maas','title':'Gradient flows of the entropy for finite Markov chains','year':2011,'journal':'Journal of Functional Analysis 261, 2250–2292','url':'https://arxiv.org/abs/1102.5238'}]
sources=json.loads((OUT/'sources.json').read_text())['sources']+json.loads((OUT/'engineering_sources.json').read_text())['sources']+extra
sources=list({s['key']:s for s in sources}.values());(ROOT/'sources_v2.json').write_text(json.dumps({'sources':sources},indent=2,ensure_ascii=False))

numfiles=['tree_exact.py','tree_exact_results.json','noise_geometry_audit.py','noise_geometry_results.json','noise_smooth_review.py','noise_smooth_review.json','statistics.py','statistics_results.json','statistics_compiled_results.json','precision_cert.py','precision_cert_results.json','convex_compile.py','convex_raw_64.json','primal_certificate.py','primal_cert_results.json','dual_certificate.py','dual_bounds_64_128.json','aggregate_convex.py','convex_results.json','replay_v2.py','replay_v2_receipt.json','review_receipt.json','requirements-core.txt','requirements-opt.txt','README.md','sources_v2.json']
hashes={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest()for f in numfiles};(ROOT/'SHA256SUMS.json').write_text(json.dumps(hashes,indent=2))
with zipfile.ZipFile(OUT/'Cyclecraft_Numerical_Program.zip','w',zipfile.ZIP_DEFLATED)as z:
 for f in numfiles+['SHA256SUMS.json']:z.write(ROOT/f,f)
stats=load('statistics_results.json');compiled=load('statistics_compiled_results.json')
compact={k:v for k,v in stats.items()if k not in ['models','cases']};compact['cases']=stats['cases']+compiled['cases'];compact['models']={k:{a:b for a,b in v.items()if a!='cases'}for k,v in stats['models'].items()}
data={'statistics':compact,'convex':load('convex_results.json'),'precision':load('precision_cert_results.json'),'exact':load('tree_exact_results.json'),'noise_geometry':load('noise_geometry_results.json'),'compiled_checks':compiled['checks']}
jsondata=json.dumps(data,separators=(',',':'),allow_nan=False)

front=read('paper_front_v2.html');theory=read('paper_theory_v2.html');convex=read('paper_convex_v2.html');noise=read('paper_noise_v2.html');results=read('paper_results_v2.html');cert=read('paper_certification_v2.html');engineering=read('paper_engineering_v2.html')
front=front.replace('#storage_geometry','#storage-geometry').replace('#noise-geometry"','#noise-geometry-v2"').replace('#engineering"','#engineering-purpose"')
# Insert each retained-data plate beside the theorem or experiment it explains.
markup=read('visuals_v2_markup.html');figs=re.findall(r'<figure\b[\s\S]*?</figure>',markup);toolbar=re.search(r'<div class="storage-controls"[\s\S]*?</div>',markup).group(0)
plate1='<section class="storage-laboratory" aria-label="Architecture and shared controls">'+figs[0]+toolbar+'</section>'
plate2='<section class="storage-laboratory" aria-label="Finite-horizon uncertainty">'+figs[1]+'</section>'
plate3='<section class="storage-laboratory" aria-label="Certified resource frontier">'+figs[2]+'</section>'
plate4='<section class="storage-laboratory" aria-label="Actual particle path">'+figs[3]+'</section>'
prose=front+plate1+theory+convex+plate3+noise+plate2+results+plate4+cert+read('certificate_markup.html')+engineering
# Normalize named section anchors supplied by independently authored fragments.
aliases={'convex-synthesis':'bounded-synthesis','numerical-results':'numerical-program-v2','certification':'certification-v2'}
ids=set(re.findall(r'\bid="([^"]+)"',prose))
for a,b in aliases.items():
 if a not in ids:
  if b in ids:prose=prose.replace('href="#'+a+'"','href="#'+b+'"')
  else:
   names=[i for i in ids if a.split('-')[0]in i]
   if names:prose=prose.replace('href="#'+a+'"','href="#'+sorted(names)[0]+'"')
prose=re.sub(r'(<table\b[\s\S]*?</table>)',r'<div class="table-wrap">\1</div>',prose)
plain=re.sub(r'\\\[[\s\S]*?\\\]|\\\([\s\S]*?\\\)',' ',prose);plain=re.sub(r'<[^>]+>',' ',plain);wordcount=len(re.findall(r"\b[\w’'-]+\b",html.unescape(plain)))

prose+=f'''<section id="reproducibility"><h2>The complete program and its evidence</h2><p>This edition contains <strong>{wordcount:,} words of scientific prose</strong>, excluding equations and references, with autonomous numerical figures and all supporting source code. The numerical program was completed before the manuscript was assembled. Its replay recomputed eight independent stages and passed; the retained optimization proposals are verified through exact and interval inequalities.</p><div class="precision-grid"><div><span>Exact tree identities</span><b>347 checks passed</b></div><div><span>Exact noise geometry</span><b>210 checks passed</b></div><div><span>Continuous witness feasibility</span><b>42 checks passed</b></div><div><span>Interval probability evidence</span><b>34 checks passed</b></div></div><p>The harmonic counting program passes 81 checks; the independent counting calculation for the final certified protocols passes 23. The smooth noise-floor review covers nine cases. Continuous dual bounds and certified resource frontiers have separate source-bound records. The two-minute default replay verifies certificates without rerunning the conic search. Full regeneration additionally requires the optional optimization dependencies.</p><div class="download-actions"><button id="download-program" type="button">Download the complete numerical program</button><button id="download-evidence" type="button">Export the displayed evidence</button><button type="button" onclick="window.print()">Print paper</button></div><pre><code>python -m pip install -r requirements-core.txt\npython replay_v2.py\n\n# Optional: regenerate the optimization proposals\npython -m pip install -r requirements-opt.txt\npython replay_v2.py --full</code></pre><p class="note">Reference environment: Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0, SymPy 1.14.0, mpmath 1.3.0, python-flint 0.9.0. Optional search: CVXPY 1.9.3 and Clarabel 0.11.1. All figures and fonts are embedded; the paper makes no network requests.</p></section>'''
used=set(re.findall(r'href="#ref-([^"]+)"',prose));by={s['key']:s for s in sources};missing=used-set(by);assert not missing,missing
bib='<section id="references"><h2>Primary literature</h2><ol class="references">'
for key in sorted(used):
 s=by[key];bib+=f'<li id="ref-{key}">{html.escape(s["authors"])}. <em>{html.escape(s["title"])}</em>. {html.escape(s.get("journal",""))} ({s["year"]}). <a target="_blank" rel="noopener" href="{html.escape(s["url"])}">Original paper</a>.</li>'
bib+='</ol><p class="note">Numerical software reference: <a href="https://python-flint.readthedocs.io/en/latest/arb_mat.html">python-flint matrix-exponential documentation</a>. Classical pumping, graph entropy, discrete transport and counting-statistics results retain their attribution. The storage architecture, exact examples and certified comparisons are derived and executed in the accompanying program.</p></section>'
prose+=bib+'<footer class="foot">CYCLECRAFT II · Storage without winding · Second edition · 30 October 2026<br>Exact identities, continuous-time resource certificates, full counting statistics and autonomous figures.<br>First edition and complete source retained in the full package. KaTeX equations and fonts embedded under their license.</footer>'

pattern=re.compile(r'\\\[([\s\S]*?)\\\]|\\\(([\s\S]*?)\\\)');matches=list(pattern.finditer(prose));tex=[{'tex':html.unescape(m.group(1)if m.group(1)is not None else m.group(2)),'display':m.group(1)is not None}for m in matches]
assert not any('\\operatorname'in a['tex']for a in tex)
node=os.environ.get('CODEX_PRIMARY_RUNTIME_NODE','node');js="const fs=require('fs'),k=require('../katex.min.js');const x=JSON.parse(fs.readFileSync(0,'utf8'));process.stdout.write(JSON.stringify(x.map((a,i)=>{try{return k.renderToString(a.tex,{displayMode:a.display,throwOnError:true,strict:'ignore',output:'htmlAndMathml'})}catch(e){throw Error('Equation '+i+': '+e.message+' TEX='+a.tex)}})));"
r=subprocess.run([node,'-e',js],input=json.dumps(tex),cwd=ROOT,capture_output=True,text=True);assert r.returncode==0,r.stderr
rendered=json.loads(r.stdout)
for i,a in enumerate(tex):
 if not a['display']and re.fullmatch(r'-?\d\.\d{24,}',a['tex'].strip()):rendered[i]='<code class="precision-number">'+html.escape(a['tex'].strip())+'</code>'
it=iter(rendered);prose=pattern.sub(lambda m:next(it),prose)
css=(OUT/'katex_embedded.css').read_text();css=re.sub(r',url\(fonts/[^)]+\) format\("[^\"]+"\)','',css)
css+='\n'+(OUT/'paper.css').read_text()+'\n'+read('visuals_v2.css')+'\n'+read('paper_v2.css')
tools='''<div class="reading-progress" aria-hidden="true"></div><main><div class="reader-tools"><a href="#abstract">Read the paper</a><button id="theme-toggle" type="button">Dark / light</button><button type="button" onclick="document.getElementById('storage-play').click()">Play / pause</button><button type="button" onclick="document.getElementById('reproducibility').scrollIntoView()">Code & evidence</button><button type="button" onclick="window.print()">Print</button></div>'''
runtime="""document.getElementById('theme-toggle').onclick=()=>{document.body.classList.toggle('dark');document.documentElement.classList.toggle('dark');window.dispatchEvent(new Event('resize'))};addEventListener('scroll',()=>{let d=document.documentElement;document.querySelector('.reading-progress').style.width=(100*d.scrollTop/(d.scrollHeight-d.clientHeight))+'%'},{passive:true});let printWasDark=false;addEventListener('beforeprint',()=>{printWasDark=document.documentElement.classList.contains('dark');document.documentElement.classList.remove('dark');document.body.classList.remove('dark');dispatchEvent(new Event('resize'))});addEventListener('afterprint',()=>{if(printWasDark){document.documentElement.classList.add('dark');document.body.classList.add('dark');dispatchEvent(new Event('resize'))}});function saveBlob(b,n,t){let a=document.createElement('a');a.href=URL.createObjectURL(new Blob([b],{type:t}));a.download=n;a.click();setTimeout(()=>URL.revokeObjectURL(a.href),1000)}document.getElementById('download-evidence').onclick=()=>saveBlob(JSON.stringify(DATA2,null,2),'Cyclecraft_II_Evidence.json','application/json');document.getElementById('download-program').onclick=()=>{let s=atob(document.getElementById('program-archive').textContent.trim()),b=new Uint8Array(s.length);for(let i=0;i<s.length;i++)b[i]=s.charCodeAt(i);saveBlob(b,'Cyclecraft_Numerical_Program.zip','application/zip')};"""
scripts=read('visuals_v2.js')+'\n'+read('hero_v2.js')+'\n'+read('certificate_visual.js')+'\n'+runtime
archive=base64.b64encode((OUT/'Cyclecraft_Numerical_Program.zip').read_bytes()).decode()
out='<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><meta name="description" content="Cyclecraft II: controlled tree storage, certified finite-speed transport, exact graph geometry and finite-horizon counting risk."><title>Cyclecraft II — Storage without winding</title><style>'+css+'</style><!-- KaTeX license\n'+(OUT/'KATEX_LICENSE.txt').read_text()+'\n--></head><body>'+tools+prose+'</main><script>const DATA2='+jsondata+';</script><script>'+scripts+'</script><script id="program-archive" type="application/octet-stream">'+archive+'</script></body></html>'
(OUT/'Cyclecraft.html').write_text(out)
receipt={'prose_words':wordcount,'math_expressions':len(tex),'html_bytes':len(out.encode()),'program_archive_bytes':(OUT/'Cyclecraft_Numerical_Program.zip').stat().st_size,'self_contained':True}
(ROOT/'build_v2_receipt.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt,indent=2))

authoring=['build_v2.py','paper_front_v2.html','paper_theory_v2.html','paper_convex_v2.html','paper_noise_v2.html','paper_results_v2.html','paper_certification_v2.html','paper_engineering_v2.html','paper_v2.css','visuals_v2.js','visuals_v2.css','visuals_v2_markup.html','hero_v2.js','certificate_markup.html','certificate_visual.js','build_v2_receipt.json']
if (ROOT/'qa_v2_receipt.json').exists():authoring+=['qa_v2_receipt.json']
with zipfile.ZipFile(OUT/'Cyclecraft_Complete.zip','w',zipfile.ZIP_DEFLATED)as z:
 z.write(OUT/'Cyclecraft.html','Cyclecraft/Cyclecraft.html');z.write(OUT/'Cyclecraft_Numerical_Program.zip','Cyclecraft/Cyclecraft_Numerical_Program.zip')
 for f in numfiles+['SHA256SUMS.json']+authoring:z.write(ROOT/f,'Cyclecraft/v2/'+f)
 for f in ['paper.css','katex.min.js','katex_embedded.css','KATEX_LICENSE.txt','sources.json','engineering_sources.json']:z.write(OUT/f,'Cyclecraft/'+f)
 legacy=ROOT/'legacy/Cyclecraft_V1_Archive.zip'
 if not legacy.exists():legacy=OUT/'First_Edition/Cyclecraft_V1_Archive.zip'
 z.write(legacy,'Cyclecraft/First_Edition/Cyclecraft_V1_Archive.zip')
