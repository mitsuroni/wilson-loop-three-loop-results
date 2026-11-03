"""Exact local Hlog coordinates for rational letters and explicitly named periods.
This extends computation-side output only; no Kernel AST or formal status changes.
"""
import ast,re
from pathlib import Path
from fractions import Fraction as Q
import sympy as sp
from check_math import ROOT,module
core=module('t4_wide_word_core',ROOT/'workflow/checks/single-opposite/hlog_algebra.py')
PERIODS=('zeta_2','zeta_3','zeta_5','log_2','A4','A6');WEIGHTS=(2,3,5,1,4,6)
ZERO=(0,)*6
core.ZERO=ZERO;core.ONE={((),ZERO):Q(1)};core.LOG={((Q(0),),ZERO):Q(1)}
def period(n):
 if n in (2,3,5):
  ps=list(ZERO);ps[(2,3,5).index(n)]=1;return {((),tuple(ps)):Q(1)}
 if n==4:return core.scale(core.power(period(2),2),Q(2,5))
 if n==6:return core.scale(core.power(period(2),3),Q(8,35))
 raise ValueError(n)
core.period=period
add,scale,mul,power,from_sympy=core.add,core.scale,core.mul,core.power,core.from_sympy

def parse(path):return parse_text(Path(path).read_text(encoding='utf-8').split('=',1)[1].strip())
def parse_text(text):
 words=[]
 def hlog(m):
  word=tuple(Q(a.strip()) for a in m[1].split(',') if a.strip());assert all(a in (0,-1,Q(-1,2),-2) for a in word),word
  if word not in words:words.append(word)
  return 'g'+str(words.index(word))
 text=re.sub(r'Hlog\(z,\[([^]]*)\]\)',hlog,text)
 text=text.replace('zeta[1,-3]','A4').replace('zeta[1, -3]','A4').replace('zeta[1,-5]','A6').replace('zeta[1, -5]','A6').replace('ln(2)','L2')
 text=re.sub(r'zeta\[(\d+)\]',r'Z\1',text).replace('^','**')
 gs=sp.symbols('g:'+str(len(words)));zs=sp.symbols('Z2 Z3 Z5 L2 A4 A6');names={str(s):s for s in (*gs,*zs)}
 names['Z4']=sp.Rational(2,5)*zs[0]**2;names['Z6']=sp.Rational(8,35)*zs[0]**3
 def node(n):
  if isinstance(n,ast.Constant) and type(n.value)is int:return sp.Integer(n.value)
  if isinstance(n,ast.Name) and n.id in names:return names[n.id]
  if isinstance(n,ast.UnaryOp) and isinstance(n.op,ast.USub):return -node(n.operand)
  if isinstance(n,ast.BinOp):
   a,b=node(n.left),node(n.right)
   if isinstance(n.op,ast.Add):return a+b
   if isinstance(n.op,ast.Sub):return a-b
   if isinstance(n.op,ast.Mult):return a*b
   if isinstance(n.op,ast.Div):return a/b
   if isinstance(n.op,ast.Pow) and b.is_Integer and b>=0:return a**b
  raise ValueError('unsupported expression node: '+ast.dump(n)[:200])
 expr=sp.expand(node(ast.parse(text,mode='eval').body));assert not gs or sp.Poly(expr,*gs).total_degree()<=1
 out={}
 for g,w in [*zip(gs,words),(None,())]:
  c=expr.subs(dict.fromkeys(gs,0)) if g is None else expr.coeff(g)
  for ps,v in sp.Poly(c,*zs).terms():
   assert v.is_Rational
   if v:out[(w,ps)]=Q(int(v.p),int(v.q))
 return out

def lift(v):return {(tuple(Q(a) for a in w),tuple(ps)+(0,0,0)):c for (w,ps),c in v.items()}
def lower(v):
 assert all(set(w)<={0,-1} and not any(ps[3:]) for w,ps in v)
 return {(tuple(int(a) for a in w),ps[:3]):c for (w,ps),c in v.items()}
def weight(key):return len(key[0])+sum(n*p for n,p in zip(WEIGHTS,key[1]))
def serialize(v):return [{'word':[int(a) if a.denominator==1 else str(a) for a in w],'period_powers':list(ps),'coefficient':str(c)} for (w,ps),c in sorted(v.items())]
def vec(rows):
 if rows and 'zeta_powers' in rows[0]:return {(tuple(Q(a) for a in r['word']),tuple(r['zeta_powers'])+(0,0,0)):Q(r['coefficient']) for r in rows}
 return {(tuple(Q(a) for a in r['word']),tuple(r['period_powers'])):Q(r['coefficient']) for r in rows}