import MDAnalysis
import MDAnalysis.analysis.rdf
import matplotlib.pyplot as plt

u = MDAnalysis.Universe('XDATCAR.pdb', permissive=True)
g1= u.select_atoms('type O')
g2= u.select_atoms('type H')
rdf = MDAnalysis.analysis.rdf.InterRDF(g1,g2,nbins=75, range=(0.0, min(u.dimensions[:3])/2.0))
           
rdf.run()

fig = plt.figure(figsize=(5,4))
ax = fig.add_subplot(111)
ax.plot(rdf.bins, rdf.rdf, 'k-',  label="rdf")

ax.legend(loc="best")
ax.set_xlabel(r"Distance ($\AA$)")
ax.set_ylabel(r"RDF")
fig.savefig("RDF_all.png")
#plt.show()