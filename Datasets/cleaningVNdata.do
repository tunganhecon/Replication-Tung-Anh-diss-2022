*this do-file replicates the data cleaning process for the raw Vietnamese data, where CPI and nominal interest rate are monthly, GDP is quarterly

clear

*cd to data folder
cd "C:\Users\Tung Anh\Documents\Redoing codes for diss\Datasets\Vietnam"

*importing CPI
import excel "C:\Users\Tung Anh\Documents\Redoing codes for diss\Datasets\Vietnam\VNCPI.xlsx", sheet("Historical Values") cellrange(A1:B324)

*cleaning CPI 
gen mdate = monthly(A, "MY")
format mdate %tm
drop A
rename B CPI
sort mdate 
gen qdate = qofd(dofm(mdate))
format qdate %tq
collapse (mean) CPI, by(qdate)
save "C:\Users\Tung Anh\Documents\Redoing codes for diss\Datasets\Vietnam\CPI.dta"

clear

*cleaning IR
import excel "C:\Users\Tung Anh\Documents\Redoing codes for diss\Datasets\Vietnam\VNIR.xlsx", sheet("Sheet1")
drop in 1/1
gen mdate = ym(A, B)
format mdate %tm
drop A B
order mdate C
rename C IR
gen qdate = qofd(dofm(mdate))
format qdate %tq
gen logIR = ln(IR)
collapse (mean) IR=logIR, by(qdate)
replace IR = exp(IR)
save "C:\Users\Tung Anh\Documents\Redoing codes for diss\Datasets\Vietnam\IR.dta"

clear

*cleaning GDP 
import excel "C:\Users\Tung Anh\Documents\Redoing codes for diss\Datasets\Vietnam\VNGDP.xlsx", sheet("Historical Values")
gen qdate = quarterly(substr(A, 2, .), "QY")
format qdate %tq
drop A
rename B GDP
sort qdate
order qdate GDP
save "C:\Users\Tung Anh\Documents\Redoing codes for diss\Datasets\Vietnam\GDP.dta"

clear
*merging datasets
use "C:\Users\Tung Anh\Documents\Redoing codes for diss\Datasets\Vietnam\GDP.dta", clear
joinby qdate using "C:\Users\Tung Anh\Documents\Redoing codes for diss\Datasets\Vietnam\CPI.dta", unmatched(both)
drop _merge
joinby qdate using "C:\Users\Tung Anh\Documents\Redoing codes for diss\Datasets\Vietnam\IR.dta", unmatched(both)
drop _merge

*converting nominal GDP to real GDP
replace GDP = GDP/CPI

erase GDP.dta
erase CPI.dta
erase IR.dta 

*generating growth rate variables
tsset qdate

tsfilter hp GDPtrend = GDP
tsfilter hp CPItrend = CPI 
gen deGDP = GDP - GDPtrend
gen deCPI = CPI - CPItrend
gen lnGDP = log(deGDP)
gen lnCPI = log(deCPI)

replace GDP =  d.lnGDP*400
replace CPI = d.lnCPI*100
gen dIR = d.IR

drop IR GDPtrend CPItrend deGDP deCPI lnGDP lnCPI 

keep if qdate > yq(2000, 4)

export excel using "VNdata", firstrow(variables) nolabel keepcellfmt