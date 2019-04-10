import MDAnalysis
import MDAnalysis.analysis.rms
import matplotlib.pyplot as plt

u = MDAnalysis.Universe('XDATCAR.pdb', permissive=True)
ref = MDAnalysis.Universe('XDATCAR.pdb', permissive=True)     # reference (with the default ref_frame=0)
ref.trajectory[0] #use first frame as reference
R = MDAnalysis.analysis.rms.RMSD(u, ref,
           select="all",         # superimpose on whole backbone of all atoms # align based on all atoms
           groupselections=["type H","type O"],
           filename="rmsd_all.dat",center=True)#,   # CORE
timestep=0.0005  #0.5fs from fs to ps as Reader has no dt information, set to 1.0 ps          
R.run()
rmsd = R.rmsd.T   # transpose makes it easier for plotting
time = rmsd[1]*timestep

fig = plt.figure(figsize=(5,4))
ax = fig.add_subplot(111)
ax.plot(time, rmsd[2], 'k-',  label="all")
ax.plot(time, rmsd[3], 'r--', label="type H")
ax.plot(time, rmsd[4], 'b--', label="type O")
ax.legend(loc="best")
ax.set_xlabel("time (ps)")
ax.set_ylabel(r"RMSD ($\AA$)")
fig.savefig("rmsd_md_analysis.png")
