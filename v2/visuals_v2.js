/* CYCLECRAFT II: retained-data numerical plates. No plotting imports. */
(() => {
 'use strict';
 const raw=typeof DATA2!=='undefined'?DATA2:(window.DATA2||{});
 const data=raw.visuals||raw;
 const stats=raw.statistics||raw.stats||raw.statistics_results||data;
 const convex=raw.convex||raw.convex_results||data;
 const ids=['storage-network','storage-noise','storage-frontier','storage-path'];
 const canvases=ids.map(id=>document.getElementById(id));
 if(canvases.some(c=>!c))return;
 const canonical=a=>a==='cyclic'?'cycle':a==='leaves'?'star':a;
 const architectureNames={cycle:'Fenced auxiliary cycle',star:'Two storage leaves',chain:'Storage chain'};
 const shortNames={cycle:'Cycle',star:'Leaves',chain:'Chain'};
 const designNames={certified_assisted:'Resource-certified source protocol · counting statistics in binary64',baseline:'Harmonic · frozen storage wells',assisted:'Assisted · common amplitude',assisted_matched:'Assisted · matched entropy'};
 let cases=Array.isArray(data.cases)?data.cases:Array.isArray(stats.cases)?stats.cases:[];
 if(!cases.length&&stats.models){Object.entries(stats.models).forEach(([architecture,model])=>Object.entries(model.cases||{}).forEach(([design,c])=>cases.push({...c,architecture,design,edges:model.edges})));}
 cases=cases.map(c=>({...c,architecture:canonical(c.architecture),mean_target:c.mean_target??c.target_mean,mean_aux:c.mean_aux??c.aux_mean,variance_target:c.variance_target??c.covariance_asymptotic?.[0]?.[0],variance_aux:c.variance_aux??c.covariance_asymptotic?.[1]?.[1],finite:c.finite||c.finite_horizons?.map(f=>({cycles:f.n_cycles,mean_target:f.mean2?.[0],mean_aux:f.mean2?.[1],variance_target:f.covariance2?.[0]?.[0],variance_aux:f.covariance2?.[1]?.[1]}))||[]}));
 let frontier=data.frontier||convex.frontier_rows||[];
 frontier=frontier.map(r=>({...r,architecture:canonical(r.architecture),budget:r.budget??r.entropy_budget,assisted_feasible:r.assisted_feasible??r.assisted_current_lower??r.assisted_current,baseline_upper:r.baseline_upper??r.baseline_continuous_upper,certified:r.certified??(Number.isFinite(r.baseline_continuous_upper)&&Number.isFinite(r.certified_gap)),strict:(r.certified_gap_lower??r.certified_gap)>0}));
 const selectDesign=document.getElementById('storage-design'),selectArchitecture=document.getElementById('storage-architecture');
 const play=document.getElementById('storage-play'),scrub=document.getElementById('storage-phase'),phaseOut=document.getElementById('storage-phase-value'),status=document.getElementById('storage-status');
 const reduced=matchMedia('(prefers-reduced-motion: reduce)');
 const state={phase:.08,design:'certified_assisted',architecture:'star',playing:!reduced.matches,last:0};
 const active=new Map(canvases.map(c=>[c,true]));
 let col;
 function colors(){const c=getComputedStyle(document.documentElement);const read=(p,d)=>c.getPropertyValue(p).trim()||d;col={ink:read('--ink','#182d38'),muted:read('--muted','#637276'),paper:read('--paper','#f7f4eb'),teal:read('--teal','#087f78'),accent:read('--accent','#b64e28'),rose:read('--rose','#a94b6a')};return col;}
 colors();
 const archColor=a=>a==='star'?col.teal:a==='chain'?col.accent:col.muted;
 function text(c,s,x,y,size=11,color=col.muted,align='left',weight=400){c.fillStyle=color;c.font=`${weight} ${size}px Arial`;c.textAlign=align;c.textBaseline='alphabetic';c.fillText(String(s),x,y);}
 function line(c,a,b,color=col.muted,width=1,alpha=1,dash=[]){c.save();c.strokeStyle=color;c.lineWidth=width;c.globalAlpha=alpha;c.setLineDash(dash);c.beginPath();c.moveTo(...a);c.lineTo(...b);c.stroke();c.restore();}
 function path(c,points,color,width=1.5,alpha=1,dash=[]){if(!points.length)return;c.save();c.strokeStyle=color;c.lineWidth=width;c.globalAlpha=alpha;c.setLineDash(dash);c.beginPath();points.forEach((p,i)=>i?c.lineTo(...p):c.moveTo(...p));c.stroke();c.restore();}
 function dot(c,x,y,r,color,alpha=1){c.save();c.fillStyle=color;c.globalAlpha=alpha;c.beginPath();c.arc(x,y,r,0,Math.PI*2);c.fill();c.restore();}
 function ring(c,x,y,r,color,width=1,dash=[]){c.save();c.strokeStyle=color;c.lineWidth=width;c.setLineDash(dash);c.beginPath();c.arc(x,y,r,0,Math.PI*2);c.stroke();c.restore();}
 function arrow(c,x,y,angle,color,size=4){c.save();c.fillStyle=color;c.translate(x,y);c.rotate(angle);c.beginPath();c.moveTo(size,0);c.lineTo(-size,-size*.6);c.lineTo(-size,size*.6);c.closePath();c.fill();c.restore();}
 function setup(canvas){const r=canvas.getBoundingClientRect(),w=Math.max(1,r.width),h=r.height,d=Math.min(devicePixelRatio||1,2.5);if(canvas.width!==Math.round(w*d)||canvas.height!==Math.round(h*d)){canvas.width=Math.round(w*d);canvas.height=Math.round(h*d);}const c=canvas.getContext('2d');c.setTransform(d,0,0,d,0,0);c.clearRect(0,0,w,h);return {c,w,h};}
 const fmt=(x,n=5)=>!Number.isFinite(x)?'—':Math.abs(x)<1e-7&&x!==0?x.toExponential(2):Math.abs(x)>=100?x.toFixed(1):x.toFixed(n);
 function currentCase(architecture=state.architecture,design=state.design){return cases.find(c=>c.architecture===architecture&&c.design===design)||cases.find(c=>c.architecture===architecture&&c.design==='assisted')||cases.find(c=>c.architecture===architecture);}
 function sample(c,phase=state.phase){const t=c?.trajectory;if(!t?.p?.length)return null;const n=t.p.length-1,z=Math.min(n-1,Math.floor(phase*n)),f=phase*n-z;const at=k=>{if(!t[k])return null;let a=t[k][z],b=t[k][z+1];return Array.isArray(a)?a.map((v,i)=>v+(b[i]-v)*f):a+(b-a)*f;};return {p:at('p'),pi:at('pi'),j:at('current'),counts:at('counts'),q:at('q')||at('conductance'),E:at('E')};}
 function unavailable(c,w,h,message='This figure requires its retained numerical dataset.'){c.font='12px Arial';let lines=[],line='';for(const word of message.split(' ')){let test=line?line+' '+word:word;if(c.measureText(test).width>w-32&&line){lines.push(line);line=word}else line=test}if(line)lines.push(line);lines.forEach((s,i)=>text(c,s,w/2,h/2+(i-(lines.length-1)/2)*20,12,col.muted,'center'));}
 function network(){const{c,w,h}=setup(canvases[0]),mobile=w<650,architectures=['cycle','star','chain'];
  if(!cases.length){unavailable(c,w,h);return;}
  const rows=mobile?3:1,cols=mobile?1:3,pw=w/cols,ph=h/rows;
  architectures.forEach((architecture,k)=>{const r=currentCase(architecture),s=sample(r);if(!r||!s)return;const ox=mobile?0:k*pw,oy=mobile?k*ph:0,px=x=>ox+x*pw,py=y=>oy+y*ph;
   text(c,shortNames[architecture].toUpperCase(),px(.03),py(.055),11,archColor(architecture),'left',600);
   if(state.architecture===architecture)line(c,[px(.03),py(.08)],[px(.21),py(.08)],archColor(architecture),2);
   text(c,'Qᴸ(T) '+fmt(r.mean_target),px(.97),py(.055),11,col.teal,'right',600);
   const positions=[[px(.52),py(.51)],[px(.14),py(.27)],[px(.14),py(.75)],[px(.86),py(.27)],[px(.86),py(.75)]];
   const currents=r.trajectory.current||[],maxj=Math.max(1e-10,...currents.flat().map(Math.abs));
   [[0,1,2],[0,3,4]].forEach((tri,i)=>{if(i&&architecture!=='cycle')return;c.save();c.globalAlpha=.035;c.fillStyle=i?col.accent:col.teal;c.beginPath();tri.forEach((v,j)=>j?c.lineTo(...positions[v]):c.moveTo(...positions[v]));c.closePath();c.fill();c.restore();});
   r.edges.forEach(([u,v],e)=>{const a=positions[u],b=positions[v],dx=b[0]-a[0],dy=b[1]-a[1],len=Math.hypot(dx,dy),ux=dx/len,uy=dy/len,r1=47*Math.sqrt(Math.max(0,s.p[u]))+4,r2=47*Math.sqrt(Math.max(0,s.p[v]))+4;const aa=[a[0]+ux*r1,a[1]+uy*r1],bb=[b[0]-ux*r2,b[1]-uy*r2],color=e<3?col.teal:col.accent,strength=Math.min(1,Math.abs(s.j?.[e]||0)/maxj),q=s.q?.[e]??1;
    line(c,aa,bb,col.muted,1+Math.min(q,3)*1.25,.15);line(c,aa,bb,color,.8+2.2*strength,.3+.6*strength);arrow(c,(aa[0]+bb[0])/2,(aa[1]+bb[1])/2,Math.atan2(dy,dx)+((s.j?.[e]||0)<0?Math.PI:0),color,3+strength);
    if(strength>.025)for(let j=0;j<2;j++){let f=(state.phase*(2+strength*3)+j/2)%1;if(s.j[e]<0)f=1-f;dot(c,aa[0]+(bb[0]-aa[0])*f,aa[1]+(bb[1]-aa[1])*f,2,color,.7);}
   });
   positions.forEach(([x,y],i)=>{let rp=47*Math.sqrt(Math.max(0,s.p[i])),req=47*Math.sqrt(Math.max(0,s.pi?.[i]??s.p[i]));dot(c,x,y,rp+4,col.paper);dot(c,x,y,rp,col.accent,.13);ring(c,x,y,req,col.muted,1,[2,3]);ring(c,x,y,rp,col.accent,1.4);text(c,i,x,y+4,12,col.ink,'center',600);const direction=i===1||i===3?-1:1;const ly=y+direction*(Math.max(rp,req)+16)+(direction>0?4:0);text(c,s.p[i].toFixed(3),x,ly,11,col.muted,'center');});
   const aux=s.counts?.[3]??0;const target=s.counts?.[0]??0;
   text(c,'count e₀ = '+fmt(target),px(.03),py(.925),11,col.teal);
   text(c,'count e₃ = '+fmt(aux),px(.97),py(.925),11,col.accent,'right');
   text(c,architecture==='cycle'?'zero mean auxiliary circulation':'bounded auxiliary signed counts',px(.5),py(.99),11,col.muted,'center');
  });
 }
 function range(v,pad=.12){let lo=Math.min(...v),hi=Math.max(...v);if(!Number.isFinite(lo)||!Number.isFinite(hi))return[0,1];if(Math.abs(hi-lo)<1e-12){lo-=.5;hi+=.5;}let d=(hi-lo)*pad;return[lo-d,hi+d];}
 function axes(c,r,xmin,xmax,ymin,ymax,xticks,yticks,xfmt=x=>fmt(x,2),yfmt=x=>fmt(x,3)){
  const X=x=>r.x+(x-xmin)/(xmax-xmin)*r.w,Y=y=>r.y+r.h-(y-ymin)/(ymax-ymin)*r.h;
  xticks.forEach(v=>{line(c,[X(v),r.y+r.h],[X(v),r.y+r.h+4],col.muted,1,.5);text(c,xfmt(v),X(v),r.y+r.h+19,11,col.muted,'center');});
  yticks.forEach(v=>{line(c,[r.x,Y(v)],[r.x+r.w,Y(v)],col.muted,.7,.15);text(c,yfmt(v),r.x-7,Y(v)+4,11,col.muted,'right');});
  line(c,[r.x,r.y+r.h],[r.x+r.w,r.y+r.h],col.muted,1,.5);return{X,Y};
 }
 function finite(c){return(c?.finite||[]).filter(f=>f.cycles>0&&Number.isFinite(f.variance_target)).sort((a,b)=>a.cycles-b.cycles);}
 function noise(){const{c,w,h}=setup(canvases[1]),mobile=w<650,list=['cycle','star','chain'].map(a=>({architecture:a,rows:finite(currentCase(a))})).filter(v=>v.rows.length);if(!list.length){unavailable(c,w,h);return;}
  const ns=[...new Set(list.flatMap(v=>v.rows.map(r=>r.cycles)))].sort((a,b)=>a-b),xmin=Math.log10(ns[0]),xmax=Math.log10(ns[ns.length-1]||ns[0])+ (ns.length===1?1:0);
  const rt=mobile?{x:55,y:46,w:w-79,h:235}:{x:60,y:46,w:w*.43-55,h:255};
  const ra=mobile?{x:55,y:405,w:w-79,h:205}:{x:w*.55+5,y:46,w:w*.45-27,h:255};
  const vals=list.flatMap(v=>v.rows.flatMap(r=>[Math.abs(r.mean_target),Math.sqrt(Math.max(r.variance_target,0))])).filter(v=>v>0);
  let lo=Math.floor(Math.log10(Math.min(...vals))),hi=Math.ceil(Math.log10(Math.max(...vals)));if(lo===hi)hi++;
  text(c,mobile?'TARGET · |MEAN| AND SD':'TARGET · MEAN AND STANDARD DEVIATION',rt.x,rt.y-24,11,col.teal,'left',600);
  const at=axes(c,rt,xmin,xmax,lo,hi,ns.map(Math.log10),Array.from({length:hi-lo+1},(_,i)=>lo+i),x=>String(Math.round(10**x)),x=>(10**x).toPrecision(1));
  const aux=list.flatMap(v=>v.rows.map(r=>Math.sqrt(Math.max(0,r.variance_aux)))).filter(v=>v>0);let al=Math.floor(Math.log10(Math.min(...aux,1))),ah=Math.ceil(Math.log10(Math.max(...aux,1)));if(al===ah)ah++;
  text(c,'AUXILIARY · STANDARD DEVIATION',ra.x,ra.y-24,11,col.accent,'left',600);
  const aa=axes(c,ra,xmin,xmax,al,ah,ns.map(Math.log10),Array.from({length:ah-al+1},(_,i)=>al+i),x=>String(Math.round(10**x)),x=>(10**x).toPrecision(1));
  line(c,[ra.x,aa.Y(0)],[ra.x+ra.w,aa.Y(0)],col.muted,1,.55,[3,4]);text(c,'bridge bound: |N| ≤ 1',ra.x+ra.w,aa.Y(0)-7,11,col.muted,'right');
  const index=Math.min(ns.length-1,Math.floor(state.phase*ns.length)),n=ns[index];
  list.forEach(({architecture,rows})=>{const color=archColor(architecture),selected=architecture===state.architecture;const points=(key,transform,ax)=>rows.filter(r=>transform(r[key])>0).map(r=>[ax.X(Math.log10(r.cycles)),ax.Y(Math.log10(transform(r[key])))]);
   path(c,points('mean_target',Math.abs,at),color,selected?2.6:1.3,selected?1:.7);path(c,points('variance_target',v=>Math.sqrt(Math.max(v,0)),at),color,selected?2:1.1,.85,[5,4]);path(c,points('variance_aux',v=>Math.sqrt(Math.max(v,0)),aa),color,selected?2.6:1.3,selected?1:.7);
   rows.forEach(r=>{dot(c,at.X(Math.log10(r.cycles)),at.Y(Math.log10(Math.abs(r.mean_target))),2.5,color);dot(c,aa.X(Math.log10(r.cycles)),aa.Y(Math.log10(Math.sqrt(Math.max(r.variance_aux,1e-30)))),2.5,color);});
  });
  [at,aa].forEach((a,i)=>{const r=i?ra:rt;line(c,[a.X(Math.log10(n)),r.y],[a.X(Math.log10(n)),r.y+r.h],col.muted,1,.35,[3,4]);text(c,'completed cycles · log scale',r.x+r.w/2,r.y+r.h+42,11,col.muted,'center');});
  const legendY=mobile?704:357;list.forEach((v,i)=>{const x=mobile?12+i*w/3:rt.x+i*125;line(c,[x,legendY],[x+18,legendY],archColor(v.architecture),2);text(c,shortNames[v.architecture],x+24,legendY+4,11,col.muted);});
  if(mobile){text(c,'Solid: |mean target count|',12,746,11,col.muted);text(c,'Dashed: target standard deviation',12,767,11,col.muted);}else text(c,'Solid: |mean target count|   ·   Dashed: target standard deviation',rt.x,388,11,col.muted);
 }
 function frontierDraw(){const{c,w,h}=setup(canvases[2]),rows=frontier.filter(r=>r.architecture===state.architecture&&Number.isFinite(r.budget)&&Number.isFinite(r.assisted_feasible)).sort((a,b)=>a.budget-b.budget);if(!rows.length){unavailable(c,w,h,'Resource-frontier calculations are not included for this architecture.');return;}
  const mobile=w<650,r={x:60,y:50,w:w-84,h:mobile?290:260},xs=rows.map(r=>r.budget),ys=rows.flatMap(r=>[r.assisted_feasible,r.baseline_upper,r.baseline_primal]).filter(Number.isFinite),[xl,xh]=range(xs,.08),[rawyl,yh]=range([...ys,0],.10),yl=Math.max(0,rawyl);
  text(c,shortNames[state.architecture].toUpperCase()+' · FIXED GATE WAVEFORM',r.x,20,11,archColor(state.architecture),'left',600);
  const ax=axes(c,r,xl,xh,yl,yh,xs.length>5?[xs[0],xs[Math.floor(xs.length/2)],xs.at(-1)]:xs,[yl,0,yh],x=>fmt(x,3),y=>fmt(y,4));
  const certified=rows.filter(r=>r.certified&&r.strict&&Number.isFinite(r.baseline_upper));
  if(certified.length>1){c.save();c.fillStyle=col.teal;c.globalAlpha=.09;c.beginPath();certified.forEach((v,i)=>i?c.lineTo(ax.X(v.budget),ax.Y(v.assisted_feasible)):c.moveTo(ax.X(v.budget),ax.Y(v.assisted_feasible)));certified.slice().reverse().forEach(v=>c.lineTo(ax.X(v.budget),ax.Y(v.baseline_upper)));c.closePath();c.fill();c.restore();}
  path(c,rows.map(v=>[ax.X(v.budget),ax.Y(v.assisted_feasible)]),col.teal,2.6);
  const upper=rows.filter(v=>Number.isFinite(v.baseline_upper));path(c,upper.map(v=>[ax.X(v.budget),ax.Y(v.baseline_upper)]),col.ink,1.7,1,[5,4]);
  const primal=rows.filter(v=>Number.isFinite(v.baseline_primal));path(c,primal.map(v=>[ax.X(v.budget),ax.Y(v.baseline_primal)]),col.muted,1.3,.7,[1,4]);
  rows.forEach(v=>dot(c,ax.X(v.budget),ax.Y(v.assisted_feasible),3,col.teal));upper.forEach(v=>ring(c,ax.X(v.budget),ax.Y(v.baseline_upper),3,col.ink));
  const idx=Math.min(rows.length-1,Math.floor(state.phase*rows.length)),s=rows[idx];line(c,[ax.X(s.budget),r.y],[ax.X(s.budget),r.y+r.h],col.muted,1,.35,[3,4]);
  text(c,'entropy budget Σmax',r.x+r.w/2,r.y+r.h+42,11,col.muted,'center');
  let y=r.y+r.h+73;line(c,[12,y-4],[32,y-4],col.teal,2.5);text(c,'Assisted feasible transport',40,y,11,col.teal);y+=21;line(c,[12,y-4],[32,y-4],col.ink,1.6,1,[5,4]);text(c,'Frozen-storage upper bound',40,y,11,col.ink);
  const gap=s.certified_gap_lower??(s.assisted_feasible-s.baseline_upper);if(s.certified&&Number.isFinite(gap)){text(c,gap>0?'Certified advantage ≥ '+(Math.floor(gap*1e7)/1e7).toFixed(7):'Bounds overlap · advantage not established',mobile?12:w-12,mobile?y+27:y,mobile?10.5:11,gap>0?col.teal:col.muted,mobile?'left':'right',600);}else{text(c,'Certificate status appears in the experiment table',12,y+27,mobile?10.5:11,col.muted);}
 }
 function component(edges,removed){const[a,b]=edges[removed],seen=new Set([b]),stack=[b];while(stack.length){const u=stack.pop();edges.forEach(([v,w],e)=>{if(e===removed)return;const n=v===u?w:w===u?v:-1;if(n>=0&&!seen.has(n)){seen.add(n);stack.push(n);}});}return{nodes:seen,isBridge:!seen.has(a)};}
 function pathDraw(){const{c,w,h}=setup(canvases[3]),r=currentCase(),sp=r?.sample_path;if(!sp){unavailable(c,w,h,'A retained seeded jump trajectory is required for this plate.');return;}
  const times=sp.time||sp.event_time,states=sp.state,counts=sp.aux_count;if(!times?.length||!states?.length||!counts?.length){unavailable(c,w,h);return;}
  const duration=sp.duration||times.at(-1),mobile=w<650,rr={x:49,y:60,w:w-72,h:130},rc={x:49,y:mobile?278:267,w:w-72,h:mobile?120:95};const X=t=>rr.x+t/duration*rr.w,Y=s=>rr.y+rr.h-s/4*rr.h;
  text(c,architectureNames[state.architecture],rr.x,18,12,archColor(state.architecture),'left',600);text(c,'seed '+sp.seed+(mobile?'':' · horizon '+fmt(duration,1)),rr.x,38,11,col.muted);
  [0,1,2,3,4].forEach(v=>{line(c,[rr.x,Y(v)],[rr.x+rr.w,Y(v)],col.muted,.7,.15);text(c,v,rr.x-9,Y(v)+4,11,col.muted,'right');});text(c,'state X(t)',rr.x,rr.y-9,11,col.muted);
  const end=state.phase*duration;let maxIndex=0;while(maxIndex+1<times.length&&times[maxIndex+1]<=end)maxIndex++;
  function step(values,ay,limit=times.length-1){let pts=[];for(let i=0;i<=limit;i++){if(i)pts.push([X(times[i]),ay(values[i-1])]);pts.push([X(times[i]),ay(values[i])]);}if(limit<times.length-1)pts.push([X(end),ay(values[limit])]);else pts.push([X(duration),ay(values[limit])]);return pts;}
  path(c,step(states,Y),col.muted,.9,.25);path(c,step(states,Y,maxIndex),col.teal,1.6);dot(c,X(end),Y(states[maxIndex]),3.5,col.teal);
  const comp=component(r.edges,3),indicator=states.map(s=>(comp.nodes.has(s)?1:0)-(comp.nodes.has(states[0])?1:0));let v=counts.slice();if(comp.isBridge)v.push(-1,1);const[lo,hi]=range(v,.18),CY=n=>rc.y+rc.h-(n-lo)/(hi-lo)*rc.h;
  text(c,'SIGNED AUXILIARY EDGE COUNT',rc.x,rc.y-20,11,col.accent,'left',600);[...new Set([Math.ceil(lo),0,Math.floor(hi)])].forEach(v=>{line(c,[rc.x,CY(v)],[rc.x+rc.w,CY(v)],col.muted,.7,.16);text(c,v,rc.x-9,CY(v)+4,11,col.muted,'right');});
  if(comp.isBridge){[-1,1].forEach(v=>line(c,[rc.x,CY(v)],[rc.x+rc.w,CY(v)],col.accent,1,.3,[3,4]));path(c,step(indicator,CY),col.ink,3,.35,[2,5]);}
  path(c,step(counts,CY),col.accent,1,.3);path(c,step(counts,CY,maxIndex),col.accent,2.2);dot(c,X(end),CY(counts[maxIndex]),3.5,col.accent);
  [rr,rc].forEach(b=>line(c,[X(end),b.y],[X(end),b.y+b.h],col.muted,.8,.35,[3,4]));
  [0,.5,1].forEach(v=>text(c,fmt(v*duration,0),X(v*duration),rc.y+rc.h+20,11,col.muted,'center'));
  text(c,'elapsed time',rc.x+rc.w/2,rc.y+rc.h+42,11,col.muted,'center');
  const err=comp.isBridge?Math.max(...counts.map((v,i)=>Math.abs(v-indicator[i]))):null;
  text(c,comp.isBridge?'count − endpoint indicator = '+fmt(err,0):(mobile?'Auxiliary winding can accumulate.':'Auxiliary cycle count may accumulate winding.'),12,h-10,11,comp.isBridge?col.teal:col.accent,'left',600);
 }
 const draws=[network,noise,frontierDraw,pathDraw];
 function render(all=false){draws.forEach((f,i)=>{if(all||active.get(canvases[i]))f();});scrub.value=state.phase;phaseOut.value=state.phase.toFixed(3);}
 function playback(){play.textContent=state.playing?'Pause motion':'Play motion';play.setAttribute('aria-pressed',String(state.playing));}
 function announce(){const c=currentCase();status.textContent=`${architectureNames[state.architecture]}. ${designNames[state.design]||state.design}. ${c?'Period '+c.period+'. ':''}${state.playing?'Autonomous playback.':'Paused for inspection.'}`;}
 play.onclick=()=>{state.playing=!state.playing;state.last=0;playback();announce();};
 selectDesign.onchange=()=>{state.design=selectDesign.value;render(true);announce();};selectArchitecture.onchange=()=>{state.architecture=selectArchitecture.value;render(true);announce();};
 scrub.oninput=()=>{state.phase=+scrub.value;state.playing=false;playback();render(true);};scrub.onchange=announce;
 reduced.addEventListener?.('change',e=>{if(e.matches){state.playing=false;playback();announce();}});
 if('IntersectionObserver'in window){const observer=new IntersectionObserver(entries=>entries.forEach(e=>{active.set(e.target,e.isIntersecting);if(e.isIntersecting)draws[canvases.indexOf(e.target)]();}),{rootMargin:'140px'});canvases.forEach(c=>observer.observe(c));}
 new ResizeObserver(()=>render(true)).observe(canvases[0].parentElement);
 new MutationObserver(()=>{colors();render(true);}).observe(document.documentElement,{attributes:true,attributeFilter:['class','data-theme','style']});
 addEventListener('beforeprint',()=>{state.playing=false;playback();render(true);});
 function tick(now){if(!state.last)state.last=now;const dt=Math.min(.06,(now-state.last)/1000);state.last=now;if(state.playing&&!document.hidden){state.phase=(state.phase+dt/22)%1;if([...active.values()].some(Boolean))render();}requestAnimationFrame(tick);}
 window.CyclecraftStorage={getState:()=>({...state}),render:()=>render(true),data:{cases,frontier},setPhase:p=>{state.phase=Math.max(0,Math.min(.999,+p));render(true);}};
 playback();render(true);announce();requestAnimationFrame(tick);
})();
