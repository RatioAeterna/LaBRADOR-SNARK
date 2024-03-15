# usage:  echo 761 4591 286 quotient | sage estimate.sage
# usage:  echo 761 4591 250 product | sage estimate.sage

# XXX UNDER: many underestimates and potential underestimates
# XXX OVER: many overestimates and potential overestimates

# ----- generic caching layer

import collections

class memoized(object):
  def __init__(self,func):
    self.func = func
    self.cache = {}
    self.__name__ = 'memoized:' + func.__name__
  def __call__(self,*args):
    if not isinstance(args,collections.Hashable):
      return self.func(*args)
    if not args in self.cache:
      self.cache[args] = self.func(*args)
    return self.cache[args]

# -----

def partialsums(x):
  sum = 0
  y = []
  for xi in x:
    sum += xi
    y += [sum]
  return y

@memoized
def enum(b):
  enum1 = 0.125*b*log(RR(b),2)-0.547*b+10.4
  enum2 = 0.1839*b*log(RR(b),2)-0.995*b+16.25
  return min(enum1,enum2)
  # XXX UNDER/OVER: misuse of asymptotics

def qenum(b):
  return 0.5*enum(b)
  # XXX UNDER/OVER: misuse of asymptotics

def sieve(b):
  return 0.29248125036*b
  # XXX UNDER/OVER: misuse of asymptotics

def sieverealcost(b):
  return 0.39624062518*b - 5
  # XXX UNDER/OVER: misuse of asymptotics

def qsieve(b):
  return 0.265*b
  # XXX UNDER/OVER: misuse of asymptotics
  # XXX UNDER: assumes instant QRAM

# XXX UNDER: 'free' options ignore cost of RAM
estimates = (
  ('nonq','sieving','free',sieve),
  ('nonq','sieving','real',sieverealcost),
  ('nonq','enumeration','free',enum),
  ('nonq','enumeration','real',enum),
  ('quantum','sieving','free',qsieve),
  ('quantum','sieving','real',sieverealcost),
  ('quantum','enumeration','free',qenum),
  ('quantum','enumeration','real',qenum),
)

@memoized
def choose(n,k):
  return RR(binomial(n,k))

@memoized
def delta(b):
  # XXX UNDER: experiments suggest delta is actually larger
  # XXX OVER: but maybe delta crosses below this for large b
  return (b*RR(pi*b)^(1/b)/RR(2*pi*e))^(1/(2*(b-1)))

# 1-(1-prob)^tries
def amplify(prob,tries):
  if tries == 1: return prob
  return -(RR(-prob).log1p()*tries).expm1()

def run(p,q,w,quotient):
  best = {}
  def setbest(tag,lgcost,notes):
    assert lgcost != NaN
    best[tag] = (lgcost,notes)

  if quotient:
    # XXX UNDER: incorrectly treats ntru prime as ntru classic
    equivalence = p
  else:
    # XXX OVER: assumes rotating t to \Z is optimal
    equivalence = 1
  # XXX OVER: considers only equivalence by rotations
  # XXX OVER: assumes independence across equivalence class

  # XXX OVER: limited force search
  lastneededforce = 0
  for force in range(0,p):
    if force > 20 and force > 2*lastneededforce: break
    if force > 20 and force % 10: continue
    if force > 100 and force % 20: continue
    if force > 400 and force % 40: continue

    probforce = choose(p-force,w)/choose(p,w)
    probforce = amplify(probforce,equivalence)

    spositions = p-force
    if not quotient: spositions += 1
  
    samples = p if quotient else 2*p

    # XXX OVER: limited m search
    for m in range(40,samples+1,10):
      d = spositions+m
  
      # ----- non-hybrid attacks

      # XXX OVER: limited scale search
      for scale in [1.0,sqrt(2*p/RR(3*w)),sqrt(2*m/RR(3*w))]:
        volume = scale^spositions*RR(q)^m
    
        # XXX OVER/UNDER: assumes average g weight
        target = sqrt(RR(w*scale^2+2*m/3))/volume^(1/d)

        # XXX OVER: limited block-size search
        for blocksize in range(40,m+1):
          # XXX OVER: experiments say smaller sizes often work
          if target <= delta(blocksize)^(2*blocksize-d)*sqrt(RR(d/blocksize)):
            lgprob = log(probforce,2)
            for quantum,eors,memcost,lgsvp in estimates:
              lgcostlattice = lgsvp(blocksize)
              lgcost = lgcostlattice - lgprob
              for hybrid in ['nonh','hybrid']:
                tag = hybrid,quantum,eors,memcost
                if not tag in best or lgcost < best[tag][0]:
                  notes = 'force %s lgprobforce %.6f m %s blocksize %s lgcostlattice %.6f delta %.6f scale %.6f' % (force,lgprob,m,blocksize,lgcostlattice,delta(blocksize),scale)
                  setbest(tag,lgcost,notes)
                  lastneededforce = force

      # XXX OVER: assumes dual attack is non-competitive

      # ----- hybrid attacks

      # XXX OVER: limited scale search
      # XXX: lengths below need adjustment if scale!=1
      scale = 1.0

      # XXX OVER: assumes that forcing does not help with hybrid
      # XXX OVER: limited m search in hybrid context
      # XXX: do not overcount equivalence if force!=0
      if m%40 == 0 and force == 0:

        # XXX: limited sigma search
        for sigma in range(0,spositions,40):

          if sigma == 0: continue
          # hybrid with no guessing falls back to non-hybrid covered above

          # XXX OVER: assumes even split is optimal
          sigma1 = ZZ(floor(sigma/2))
          sigma2 = ZZ(ceil(sigma/2))

          cost0 = [2^i*choose(sigma,i) for i in range(sigma+1)]
          cost1 = [2^i*choose(sigma1,i) for i in range(sigma1+1)]
          cost2 = [2^i*choose(sigma2,i) for i in range(sigma2+1)]
          prob0 = [choose(sigma,i)*choose(p-sigma,w-i)/choose(p,w) for i in range(sigma+1)]

          qcost0 = [2^i*choose(sigma,i)*(choose(p-sigma,w-i)/(2^i*choose(p,w)))^(2/3) for i in range(sigma+1)]
          # 2^i*choose(sigma,i) strings s of weight i
          # p_s = choose(p-sigma,w-i)/(2^i*choose(p,w)*probsearch)
          # sanity check: sum_i 2^i*choose(sigma,i)*p_s = 1
          # cost of search is (sum p_s^(2/3))^(3/2)
          # which is 3/2 power of:
          #  sum_i 2^i*choose(sigma,i)*choose(p-sigma,w-i)^(2/3)/(2^i*choose(p,w)*probsearch)^(2/3)
          # which is (sum_i qcost0[i])^(3/2)/probsearch

          cost0sum = partialsums(cost0)
          cost1sum = partialsums(cost1)
          cost2sum = partialsums(cost2)
          prob0sum = partialsums(prob0)
          qcost0sum = partialsums(qcost0)

          prob12 = {(i1,i2):
                    choose(sigma1,i1)*choose(sigma2,i2)*choose(p-sigma1-sigma2,w-i1-i2)/choose(p,w)
                   for i1 in range(sigma1+1)
                   for i2 in range(sigma2+1)}
          probL = [sum(prob12[i1,j] for i1 in range(j+1))
                   + sum(prob12[j,i2] for i2 in range(j))
                   for j in range(sigma1+1)]
          probLsum = partialsums(probL)

          minid = d-sigma
          b = minid-m
          T = RealDistribution('beta',((minid-1)/2,1/2))

          # XXX OVER: limited blocksize search
          for blocksize in range(40,minid,40):
            k = min(minid,floor(sqrt(RR(b*log(q))/log(delta(blocksize)))))
            lengths = [RR(q)]*(minid-k)
            next = RR(q)^(1-b/k)*delta(blocksize)^(k-1)
            for loop in range(k):
              lengths += [next]
              next /= delta(blocksize)^2

            veclen = sqrt(RR(b*w/p+2*m/3))
            # XXX UNDER/OVER: takes average weights
            # XXX UNDER/OVER: ignores anti-correlation with searched weight

            probnp = prod(1-T.cum_distribution_function(1-(g/(2*veclen))^2) for g in lengths)
            # XXX UNDER/OVER: need more experimental evidence

            if probnp == 0: continue
            probnp = amplify(probnp,equivalence)

            # option 1: simple search
            # XXX OVER: limited imax search
            for imax in range(0,sigma+1,5):
              probsearch = prob0sum[imax]
              if probsearch == 0: continue

              costsearch = cost0sum[imax]
              # XXX UNDER: ignores cost of inner loop

              for quantum,eors,memcost,lgsvp in estimates:
                lgcostlattice = lgsvp(blocksize)
                costlattice = 2^lgcostlattice
                prob = probnp * probsearch
                lgcost = log(costlattice + costsearch,2) - log(prob,2)
                tag = 'hybrid',quantum,eors,memcost
                if not tag in best or lgcost < best[tag][0]:
                  notes1 = 'm %s blocksize %s lgcostlattice %.6f delta %.6f sigma %s lgprobnp %.6f' % (
                            m,blocksize,lgcostlattice,delta(blocksize),sigma,log(probnp,2)
                           )
                  notes = notes1 + ' simple imax %s lgprobsearch %.6f lgcostsearch %.6f' % (imax,log(probsearch,2),log(costsearch,2))
                  setbest(tag,lgcost,notes)

              if probsearch > 0.99: break

            # option 2: quantum search
            # XXX OVER: limited imax search
            for imax in range(0,sigma+1,5):
              probsearch = prob0sum[imax]
              if probsearch == 0: continue

              costsearch = qcost0sum[imax]^(3/2)/probsearch
              # XXX UNDER: ignores cost of inner loop

              for quantum,eors,memcost,lgsvp in estimates:
                if quantum != 'quantum': continue
                lgcostlattice = lgsvp(blocksize)
                costlattice = 2^lgcostlattice
                prob = probnp * probsearch
                lgcost = log(costlattice + costsearch,2) - log(prob,2)
                tag = 'hybrid',quantum,eors,memcost
                if not tag in best or lgcost < best[tag][0]:
                  notes1 = 'm %s blocksize %s lgcostlattice %.6f delta %.6f sigma %s lgprobnp %.6f' % (
                            m,blocksize,lgcostlattice,delta(blocksize),sigma,log(probnp,2)
                           )
                  notes = notes1 + ' qsearch imax %s lgprobsearch %.6f lgcostsearch %.6f' % (imax,log(probsearch,2),log(costsearch,2))
                  setbest(tag,lgcost,notes)

              if probsearch > 0.99: break

            # option 3: meet-in-the-middle
            # XXX OVER: limited imax search
            for imax in range(0,sigma1+1,5):
              probmitm = probLsum[imax]
              if probmitm == 0: continue

              costmitm = cost1sum[imax] + cost2sum[imax]
              # XXX UNDER: ignores cost of inner loop

              # XXX UNDER: ignores collision probability
              
              for quantum,eors,memcost,lgsvp in estimates:
                lgcostlattice = lgsvp(blocksize)
                costlattice = 2^lgcostlattice
                prob = probnp * probmitm
                if memcost == 'real':
                  lgcost = log(costlattice + RR(costmitm)^(3/2),2) - 5
                else:
                  lgcost = log(costlattice + costmitm,2)
                lgcost -= log(prob,2)
                tag = 'hybrid',quantum,eors,memcost
                if not tag in best or lgcost < best[tag][0]:
                  notes1 = 'm %s blocksize %s lgcostlattice %.6f delta %.6f sigma %s lgprobnp %.6f' % (
                            m,blocksize,lgcostlattice,delta(blocksize),sigma,log(probnp,2)
                           )
                  notes = notes1 + ' mitm imax %s lgprobmitm %.6f lgcostmitm %.6f' % (imax,log(probmitm,2),log(costmitm,2))
                  setbest(tag,lgcost,notes)

              if probmitm > 0.99: break
  
  for quantum,eors,memcost,lgsvp in estimates:
    for hybrid in ['nonh','hybrid']:
      lgcost,notes = best[hybrid,quantum,eors,memcost]
      qorp = 'quotient' if quotient else 'product'
      sys.stdout.write('q %s p %s w %s %s %s %s %s %s lgcost %s %s\n'
        % (q,p,w,qorp,hybrid,quantum,eors,memcost,lgcost,notes))
  sys.stdout.flush()

for line in sys.stdin:
  line = line.split()
  p,q,w,qorp = ZZ(line[0]),ZZ(line[1]),ZZ(line[2]),line[3]
  assert qorp in ['product','quotient']
  if not p.is_prime():
    sys.stdout.write('warning: ntru prime requires p to be prime\n')
  if not q.is_prime():
    sys.stdout.write('warning: ntru prime requires q to be prime\n')
  Fq.<x> = GF(q)[]
  if not (x^p-x-1).is_irreducible():
    sys.stdout.write('warning: ntru prime requires x^p-x-1 irreducible mod q\n')
  # XXX: limit range of w?
  run(p,q,w,qorp == 'quotient')
  # XXX: could do both product and quotient for about the same cost
