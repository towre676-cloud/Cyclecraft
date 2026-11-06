/* A graph embedding whose height is the actual well energy, not a device drawing. */
(() => {'use strict';const can=document.getElementById('storage-hero');if(!can)return;
const cases=DATA2.statistics.cases;const states=()=>window.CyclecraftStorage?.getState()||{phase:.08,architecture:'star',design:'assisted_matched',playing:true};
const reduced=matchMedia('(prefers-reduced-motion:reduce)').matches;let visible=true,turn=.3,last=0;
new IntersectionObserver(es=>visible=es[0].isIntersecting,{rootMargin:'200px'}).observe(can);
function draw(){const st=states(),r=cases.find(c=>c.architecture===st.architecture&&c.design===st.design)||cases[0];if(!r)return;const t=r.trajectory,n=t.p.length-1,z=Math.min(n-1,Math.floor(st.phase*n)),f=st.phase*n-z;const a=k=>t[k][z].map((v,i)=>v+f*(t[k][z+1][i]-v));const p=a('p'),E=a('E'),j=a('current'),q=a('q');
const box=can.getBoundingClientRect(),w=box.width,h=box.height,d=Math.min(devicePixelRatio||1,2);if(can.width!==Math.round(w*d)||can.height!==Math.round(h*d)){can.width=Math.round(w*d);can.height=Math.round(h*d)}const c=can.getContext('2d');c.setTransform(d,0,0,d,0,0);c.fillStyle='#102831';c.fillRect(0,0,w,h);
const base=[[0,0],[-1.25,-.8],[-1.25,.8],[1.25,-.8],[1.25,.8]],scale=Math.min(w*.25,h*.31),theta=turn;
const pos=base.map(([x,z],i)=>{let u=x*Math.cos(theta)+z*Math.sin(theta),v=-x*Math.sin(theta)+z*Math.cos(theta);return[w*.5+u*scale,h*.51+v*.52*scale-E[i]*scale*.3,v]});
const text=(s,x,y,color='#c9d7d4',sz=11,al='left')=>{c.fillStyle=color;c.font=sz+'px Arial';c.textAlign=al;c.fillText(s,x,y)};
text('STORAGE WITHOUT WINDING',14,24,'#b7cdca',w<500?10:12);text(st.architecture==='star'?'TWO STORAGE LEAVES':st.architecture==='chain'?'STORAGE CHAIN':'FENCED AUXILIARY CYCLE',w-14,43,'#72cdb9',w<500?9:11,'right');
r.edges.forEach(([i,k],e)=>{let A=pos[i],B=pos[k],col=e<3?'#67c9b2':'#ca9b72';c.strokeStyle=col;c.globalAlpha=.3+Math.min(.5,Math.abs(j[e])*8);c.lineWidth=1.5+Math.min(4,q[e]);c.beginPath();c.moveTo(A[0],A[1]);c.lineTo(B[0],B[1]);c.stroke();c.globalAlpha=1;let u=(st.phase*3+(e*.13))%1;if(j[e]<0)u=1-u;c.fillStyle=col;c.beginPath();c.arc(A[0]+(B[0]-A[0])*u,A[1]+(B[1]-A[1])*u,2.5,0,Math.PI*2);c.fill()});
[0,1,2,3,4].sort((i,k)=>pos[i][2]-pos[k][2]).forEach(i=>{const[x,y]=pos[i],radius=(w<500?37:48)*Math.sqrt(p[i]);let g=c.createRadialGradient(x-radius*.3,y-radius*.4,1,x,y,radius);g.addColorStop(0,i<3?'#a8dfce':'#e2c2a2');g.addColorStop(1,i<3?'#1f6a67':'#725645');c.fillStyle=g;c.beginPath();c.arc(x,y,radius,0,Math.PI*2);c.fill();c.strokeStyle='#d6e5d2';c.lineWidth=.6;c.stroke();text(String(i),x,y+4,'#f0f5e8',12,'center');text('p='+p[i].toFixed(3),x,y+radius+17,'#a9c3be',w<500?9:11,'center')});
text('Height: Eᵢ / kBT · node area: occupation',14,h-16,'#9bb8b6',w<500?9:11);text('phase '+st.phase.toFixed(3),w-14,h-16,'#9bb8b6',10,'right')}
function frame(t){if(t-last>80){if(visible){if(states().playing&&!reduced)turn=.3+.12*Math.sin(t/14000);draw()}last=t}requestAnimationFrame(frame)}requestAnimationFrame(frame);window.addEventListener('resize',draw);draw();
})();
