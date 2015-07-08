import os


if os.path.exists("yeast_annotations_all.pkl") == False:
    import msdas
    from msdas import *
    m = AnnotationsYeast(get_yeast_filenames("all"), "yeast")   

    m.set_annotations()
    m.annotations.to_pickle("./yeast_annotations_all.pkl")


