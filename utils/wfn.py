import shutil
from interface.prepare import wfnfile, wfn_suffix, flag_wfn, soft

__all__ = ["save_wavefunction", "save_wfn"]


def save_wavefunction(time: int):
    """
    save q1 q2 q3 wave-functions
    """
    c = wfn_suffix[soft]
    if time == 0:
        shutil.copy(wfnfile, "q1." + c)
    elif time == 1:
        shutil.copy(wfnfile, "q2." + c)
    elif time == 2:
        shutil.copy(wfnfile, "q3." + c)
    else:
        shutil.copy("q2." + c, "q1." + c)
        shutil.copy("q3." + c, "q2." + c)
        shutil.copy(wfnfile, "q3." + c)


def save_wfn(time: float()):
    """
    save every nloop wave-functions
    """
    if flag_wfn:
        new_wfn = soft + "-%.2f." % float(time) + wfn_suffix[soft]
        shutil.copy(wfnfile, new_wfn)