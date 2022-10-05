## 2022-10-05
- minor bugfixes regarding Worklist.save and Sample.fit

## 2022-10-04
- fixed bug introduced by the use of pathlib.Path
- cleaned up settings.ini
- fixed encoding error

## 2022-09-28
- improved download of supplementary data

## 2022-09-26
- added changelog
- added Sample.raw to access raw-data
- removed chemical formulas from settings -> use of chempy.Substance.latex_name
- added Sample.baseline to access baseline-data
- deleted old IO-functions
- updated paths to use Path from pathlib