



m.save2midas("MD-simple-raw.csv")
m.scale_min_max_across_experiments()
m.save2midas("MD-simple-min-max_exp.csv")

e.df = e.df[m.df.columns]
e.df.fillna(0.5, inplace=True)
e.save2midas("MD-simple-errors.csv")




