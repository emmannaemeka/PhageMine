# Bundled PHANOTATE backend

PhageMine Desktop release-candidate installers include **PHANOTATE 1.6.7** as
a separate executable process for core gene prediction. PHANOTATE was written
by Katelyn McNair and contributors and is available from
<https://github.com/deprekate/PHANOTATE>.

PHANOTATE, `fastpath`, and `genbank` are distributed under the GNU General
Public License version 3 or later. Their complete license files are copied from
the installed upstream distributions into the desktop application's
`licenses/` directory at build time. The corresponding, unmodified source
distributions are available from their linked upstream projects and PyPI.

PhageMine invokes this backend as a separate process. No PHANOTATE algorithm,
default, or scientific threshold is changed. PhageMine itself remains under
its repository MIT license; the installer is an aggregate containing works
under their respective licenses.
