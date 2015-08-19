def EstimateS(SiteList):
    m = len(SiteList)
    m_inf = 0
    SpDict = {}

    for site in SiteList:
        if min(site) <= 10: m_inf += 1

        for sp in site:
            if sp in SpDict:
                SpDict[sp] += 1
            else: SpDict[sp] = 1

    IncVals = SpDict.values()
    S = len(IncVals)
    qs = [0]*10

    for i, q in enumerate(qs):
        qs[i] = IncVals.count(i+1)

    # Chao2
    q1 = qs[0]
    q2 = qs[1]
    chao2 = S + (((m-1)/m) * ((q1*(q1-1)) / (2*(q2+1))))

    var = 'und'
    if q1 > 0 and q2 > 0:
        var = q2 * (0.5*(q1/q2)**2 + (q1/q2)**3 + 0.25*(q1/q2)**4)

    # ICE
    num = 0
    n_inf = 0
    for i, qk in enumerate(qs):
        num += (i+1)*i*qk
        n_inf += (i+1)*qk

    ci = 1 - (q1/n_inf)
    gamma = (sum(qs)/ci) * (m_inf/(m_inf-1)) * (num/(n_inf**2)) - 1
    cv = max(0, gamma)
    ice = (S-sum(qs)) + (sum(qs)/ci) + ((q1/ci) * cv)

    return [chao2, var, ice, S]




""" Plot 3 """
"""
fig.add_subplot(2, 3, 3)

RADs = []
OrC = 'open'

path = mydir2 +'data/micro/EMP'+OrC
IN = path + '/EMP' + OrC + '-SbyS.txt'
num_lines = sum(1 for line in open(IN))
print 'number of sites:', num_lines

SampSizes = [16, 24, 32, 48]# , 64, 96, 128, 182, 256, 384, 512, 768, 1024, 1536, 2048, 3072, 4096, 6142, 8192, 12288, num_lines]

AvgChao = []
AvgICE = []
AvgS = []

ciLowChao = []
ciHiChao = []
ciLowICE = []
ciHiICE = []
ciLowS = []
ciHiS = []

for n in SampSizes:

    Chao2s = []
    ICEs = []
    Ss = []

    for i in range(2):

        SiteData = []
        lines = random.sample(range(1, num_lines+1), n)
        for i, line in enumerate(lines):

            data = linecache.getline(IN, line)
            data = eval(data)

            slist = data[1]
            SiteData.append(slist)

        chao2, var, ice, S = EstimateS(SiteData)

        Chao2s.append(chao2)
        ICEs.append(ice)
        Ss.append(S)

    AvgChao.append(np.mean(Chao2s))
    AvgICE.append(np.mean(ICEs))
    AvgS.append(np.mean(Ss))

    ciLowChao.append(np.percentile(Chao2s, 2.5))
    ciHiChao.append(np.percentile(Chao2s, 97.5))
    ciLowICE.append(np.percentile(ICEs, 2.5))
    ciHiICE.append(np.percentile(ICEs, 97.5))
    ciLowS.append(np.percentile(Ss, 2.5))
    ciHiS.append(np.percentile(Ss, 97.5))

    print'n:',n, ' obsS:', AvgS[-1], ' Avg Chao2:',AvgChao[-1], ' Avg ICE:',AvgICE[-1]


fig = plt.figure()
fig.add_subplot(1, 1, 1)

plt.plot(SampSizes, AvgS, color='grey', ls='-', lw=2, label='obs S')
plt.plot(SampSizes, AvgChao, color='c', ls='-', lw=2, label='Chao')
plt.plot(SampSizes, AvgICE, color='m', ls='-', lw=2, label='ICE')

plt.fill_between(SampSizes, ciLowS, ciHiS, color='grey', alpha=0.3)
plt.fill_between(SampSizes, ciLowChao, ciHiChao, color='c', alpha=0.3)
plt.fill_between(SampSizes, ciLowICE, ciHiICE, color='m', alpha=0.3)

plt.xlabel('Number of samples', fontsize=16)
plt.ylabel('Number of species', fontsize=16)
leg = plt.legend(loc=4,prop={'size':14})
leg.draw_frame(False)
"""


"""
ols = smf.ols('S ~ N', d).fit() # linear regression on log10(S) and log10(N)
ols_ci = ols.conf_int().ix['N'].tolist()
olsd = dict(a = ols.params['Intercept'], b = ols.params['N'], lb = ols_ci[0], ub = ols_ci[1])

print ols.summary()
x = np.arange(float(d.N.min()), float(d.N.max()), 0.1)
get_y = lambda a, b: a + b * x

y = get_y(olsd['a'], olsd['b'])
plt.plot(x, y, color='red', label='OLS')

plt.scatter(d.N, d.S, color='b', alpha=0.3)
plt.xlabel('log10(N)', fontsize=16)
plt.ylabel('log10(S)', fontsize=16)
"""
