# dhva_skeaf_data_analysis
This is software which can analyze the results from skeaf, a dhva calculating open-source programme. By analyzing the results files, the software can draw de orbit outline, extreme value, and the orbit in the fermi-surface. The software also achieve almost all fermi-surface painting functions.

## Environment setting
Before running python demo.py, you should get your environment ready, i work this software at my computer as settings below:
  1. Python 3.8.13
  2. PyQt5  5.15.4
  3. traits 6.3.2
  4. traitsui 7.3.1
  5. pyface 7.4.1
  6. mayavi 4.7.4
  7. numpy  1.23.0
  8. pandas 1.4.3
  9. scipy  1.8.1
  10. scikit-image  0.19.3
  11. hdbscan 0.8.28
  12. matplotlib  3.2.2

My way to build such environment is to use conda:

`conda create -n skeaf python=3.8`

`conda activate skeaf`

`conda install mayavi`

`conda install pandas`

`conda install scikit-image`

`conda install hdbscan`

`conda install matplotlib`

After doing the above creating and installing steps, your environment should be set well.


## Running the software
It is pretty easy to run the software if you have get your environment ready, before running the code you should put all files in the same directory, run it in the terminal,

`python demo.py`


## Introductions of basic functions
### Basic fermi-surface painting functons
The software achieved almost all basic fuctions of fermi-surface plotting functions. After you open the software successfully, you will see the window below,

![1.png](https://github.com/wentmil0705/dhva_skeaf_data_analysis/tree/main/test/1.png)

Then you can click **File -> Import bxsf file** to input the fermi-surface file, the basic format of bxsf file is showed in the test file, remember to keep the same format or you will get different wrong reports and the software will quit. Choose different ways to show the fermi-surface,
  1. Mode: basically simple if you only input bxsf file
  2. Interpol ratio: how good quality the picture show
  3. Brilliou zone: show in first brilliou zone or primitive brilliou zone
  4. Show slice: if you set the tip coordinate under reciprocal vetors of the vector in **section-v**, it will show the horizontal slice of the vector in fermi-surface
  5. BZ number: expand the cell
  6. Rotate: **it doesn't work now**
  7. Line: check it if you want to show the cell line
  8. Inner-outer: check it if you want to differ the inner and outer side color
  9. Axes: check it if you want to show axes
  10. More details can be set through **...**
  11. Remember to push **Updata** button after you set the parameters

Below is an example of way to show fermi surface,

![2.png](https://github.com/wentmil0705/dhva_skeaf_data_analysis/tree/main/test/2.png)

### Showing skeaf results functions
If you want to analyze the results generate from skeaf, you can click **File->Import all results file** and choose **bxsf file, results_long.out, results_orbitoutlines_invAng.out/results_orbitoutlines_invau.out**, then clik **Import**, the lower part is a table showing the calculating results, you can click one of them to show the extreme value and orbit outline, you can also check the trait checkbox and choose more than one frequencies to show the orbit in the fermi surface,
  1. Full trait: check it if you want to show the full orbits or show the part in cell
  2. Full trait full fs: check it if you want to show the full orbit in fermi-surface(by expanding cell automaticly)
  3. Mode: Seperate mode to show only the consistent part of fermi-surface with orbit on it
  4. Remember to push **Updata** button after you set the parameters

Below are some examples,

![3.png](https://github.com/wentmil0705/dhva_skeaf_data_analysis/tree/main/test/3,png)
![4.png](https://github.com/wentmil0705/dhva_skeaf_data_analysis/tree/main/test/4.png)
![5.png](https://github.com/wentmil0705/dhva_skeaf_data_analysis/tree/main/test/5.png)
![6.png](https://github.com/wentmil0705/dhva_skeaf_data_analysis/tree/main/test/6.png)
![7.png](https://github.com/wentmil0705/dhva_skeaf_data_analysis/tree/main/test/7.png)
![8.npg](https://github.com/wentmil0705/dhva_skeaf_data_analysis/tree/main/test/8.png)


## Contact me
If you get any problems, please feel free to contact me,

**MG21340029@smail.nju.edu.cn**

**532239580@qq.com**


