/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
--------------------------------------------------------------------*/

#include"lexer.h"

void lexer::ctrlrecv()
{
    int n;
	
	int ii,dd;
    
    ii=dd=0;
    
    A10 = ictrl[ii];
	ii++;
	A209 = ictrl[ii];
	ii++;
    A210 = ictrl[ii];
	ii++;
	A211 = ictrl[ii];
	ii++;
	A212 = ictrl[ii];
	ii++;
    A214 = ictrl[ii];
	ii++;
    A215 = ictrl[ii];
	ii++;
    A216 = ictrl[ii];
	ii++;
    A219 = ictrl[ii];
	ii++;
	A220 = ictrl[ii];
	ii++;
	A221 = ictrl[ii];
	ii++;
    A222 = ictrl[ii];
	ii++;
    A223 = dctrl[dd];
	dd++;
	A230 = ictrl[ii];
	ii++;
    A240 = ictrl[ii];
	ii++;
	A241 = ictrl[ii];
	ii++;
	A242 = ictrl[ii];
	ii++;
    A243 = ictrl[ii];
	ii++;
    A244 = ictrl[ii];
	ii++;
    A244_val = dctrl[dd];
	dd++;
    A245 = ictrl[ii];
	ii++;
    A245_val = dctrl[dd];
	dd++;
    A246 = ictrl[ii];
	ii++;
    A247 = dctrl[dd];
	dd++;
    A248 = ictrl[ii];
	ii++;
    A249 = dctrl[dd];
	dd++;
    A260 = ictrl[ii];
	ii++;
    A300 = ictrl[ii];
	ii++;
    A310 = ictrl[ii];
	ii++;
    A311 = ictrl[ii];
	ii++;
    A312 = ictrl[ii];
	ii++;
    A313 = ictrl[ii];
	ii++;
    A320 = ictrl[ii];
	ii++;
    A321 = ictrl[ii];
	ii++;
    A410 = ictrl[ii];
	ii++;
    A440 = dctrl[dd];
	dd++;

    M10 = ictrl[ii];
	ii++;

    I10 = ictrl[ii];
	ii++;
    I11 = ictrl[ii];
	ii++;
    I12 = ictrl[ii];
	ii++;
    I13 = ictrl[ii];
	ii++;
    I20 = ictrl[ii];
	ii++;
    I21 = ictrl[ii];
	ii++;
	I30 = ictrl[ii];
	ii++;
	I40 = ictrl[ii];
	ii++;
	I41 = ictrl[ii];
	ii++;
    I55 = dctrl[dd];
	dd++;
    I56 = ictrl[ii];
	ii++;
    I58_1 = dctrl[dd];
	dd++;
    I58_2 = dctrl[dd];
	dd++;
    I230 = ictrl[ii];
	ii++;
    I231 = dctrl[dd];
	dd++;
    I232 = dctrl[dd];
	dd++;
    I240 = ictrl[ii];
	ii++;
    I241 = dctrl[dd];
	dd++;
    I242 = dctrl[dd];
	dd++;
	
	
    B10 = ictrl[ii];
	ii++;
    B19 = ictrl[ii];
	ii++;
    B20 = ictrl[ii];
	ii++;
	B26 = ictrl[ii];
	ii++;
	B28 = ictrl[ii];
	ii++;
	B29 = dctrl[dd];
	dd++;
    B30 = ictrl[ii];
	ii++;
    B50 = dctrl[dd];
	dd++;
	B51 = dctrl[dd];
	dd++;
	B52 = dctrl[dd];
	dd++;
	B53 = dctrl[dd];
	dd++;
	B54 = dctrl[dd];
	dd++;
	B55 = dctrl[dd];
	dd++;
	B56 = dctrl[dd];
	dd++;
    B60 = ictrl[ii];
	ii++;
    B61 = ictrl[ii];
	ii++;
    B62 = ictrl[ii];
	ii++;
    B63 = ictrl[ii];
	ii++;
    B64 = ictrl[ii];
	ii++;
    B65 = dctrl[dd];
	dd++;
	B66_1 = dctrl[dd];
	dd++;
	B66_2 = dctrl[dd];
	dd++;
    B67 = ictrl[ii];
	ii++;
	B68 = ictrl[ii];
	ii++;
	B69 = ictrl[ii];
	ii++;
	B70 = ictrl[ii];
	ii++;
	B71 = ictrl[ii];
	ii++;
	B74 = ictrl[ii];
	ii++;
    B75 = ictrl[ii];
	ii++;
	B76 = ictrl[ii];
	ii++;
	B77 = ictrl[ii];
	ii++;
	B78 = ictrl[ii];
	ii++;
	B79 = dctrl[dd];
	dd++;
	B80 = dctrl[dd];
	dd++;
	B81 = ictrl[ii];
	ii++;
	B81_1 = dctrl[dd];
	dd++;
	B81_2 = dctrl[dd];
	dd++;
    B82 = ictrl[ii];
	ii++;
	B83 = dctrl[dd];
	dd++;
	B84 = ictrl[ii];
	ii++;
    B85 = ictrl[ii];
	ii++;
	B86 = ictrl[ii];
	ii++;
	B87 = ictrl[ii];
	ii++;
	B87_1 = dctrl[dd];
	dd++;
	B87_2 = dctrl[dd];
	dd++;
    B88 = dctrl[dd];
	dd++;
	B89 = ictrl[ii];
	ii++;
	B90 = ictrl[ii];
	ii++;
    B91 = ictrl[ii];
	ii++;
    B91_1 = dctrl[dd];
	dd++;
    B91_2 = dctrl[dd];
	dd++;
    B91_3 = dctrl[dd];
	dd++;
    B92 = ictrl[ii];
	ii++;
    B93 = ictrl[ii];
	ii++;
    B93_1 = dctrl[dd];
	dd++;
    B93_2 = dctrl[dd];
	dd++;
    B93_3 = dctrl[dd];
	dd++;
    B96_1 = dctrl[dd];
	dd++;
    B96_2 = dctrl[dd];
	dd++;
    B96_3 = dctrl[dd];
	dd++;
    B97 = dctrl[dd];
	dd++;
    B98 = ictrl[ii];
	ii++;
    B99 = ictrl[ii];
	ii++;
    B101 = ictrl[ii];
	ii++;
    B102 = dctrl[dd];
	dd++;
	B103 = dctrl[dd];
	dd++;
	B104 = dctrl[dd];
	dd++;
	B105_1 = dctrl[dd];
	dd++;
	B105_2 = dctrl[dd];
	dd++;
	B105_3 = dctrl[dd];
	dd++;
	B106 = ictrl[ii];
	ii++;
	B107 = ictrl[ii];
	ii++;
    B108 = ictrl[ii];
	ii++;
	B109 = ictrl[ii];
	ii++;
	B110_d = dctrl[dd];
	dd++;
	B110 = ictrl[ii];
	ii++;
	B111_zs = dctrl[dd];
	dd++;
	B111_ze = dctrl[dd];
	dd++;
    B112_zs = dctrl[dd];
	dd++;
    B112_z2 = dctrl[dd];
	dd++;
    B112_ze = dctrl[dd];
	dd++;
    B115 = ictrl[ii];
	ii++;
    B116 = ictrl[ii];
	ii++;
    B117 = dctrl[dd];
	dd++;
	B118 = dctrl[dd];
	dd++;
    B119 = dctrl[dd];
	dd++;
    B120 = dctrl[dd];
	dd++;
	B121 = ictrl[ii];
	ii++;
	B122 = dctrl[dd];
	dd++;
    B123 = dctrl[dd];
	dd++;
    B126 = dctrl[dd];
	dd++;
    B130 = ictrl[ii];
	ii++;
    B131 = dctrl[dd];
	dd++;
    B132_s = dctrl[dd];
	dd++;
    B132_e = dctrl[dd];
	dd++;
    B133 = ictrl[ii];
	ii++;
    B136 = ictrl[ii];
	ii++;
    B134 = dctrl[dd];
	dd++;
    B135 = dctrl[dd];
	dd++;
    B140_1 = dctrl[dd];
	dd++;
    B140_2 = dctrl[dd];
	dd++;
    B140_3 = dctrl[dd];
	dd++;
    B160 = ictrl[ii];
	ii++;
    B180 = ictrl[ii];
	ii++;
    B181_1 = dctrl[dd];
	dd++;
    B181_2 = dctrl[dd];
	dd++;
    B181_3 = dctrl[dd];
	dd++;
    B182_1 = dctrl[dd];
	dd++;
    B182_2 = dctrl[dd];
	dd++;
    B182_3 = dctrl[dd];
	dd++;
    B183_1 = dctrl[dd];
	dd++;
    B183_2 = dctrl[dd];
	dd++;
    B183_3 = dctrl[dd];
	dd++;
	B191 = ictrl[ii];
	ii++;
	B191_1 = dctrl[dd];
	dd++;
	B191_2 = dctrl[dd];
	dd++;
	B191_3 = dctrl[dd];
	dd++;
	B191_4 = dctrl[dd];
	dd++;
	B192 = ictrl[ii];
	ii++;
	B192_1 = dctrl[dd];
	dd++;
	B192_2 = dctrl[dd];
	dd++;
	B192_3 = dctrl[dd];
	dd++;
	B192_4 = dctrl[dd];
	dd++;
	B194_s = dctrl[dd];
	dd++;
	B194_e = dctrl[dd];
	dd++;
    B210 = ictrl[ii];
	ii++;
	B240 = ictrl[ii];
	ii++;
	B241 = ictrl[ii];
	ii++;
	B242 = ictrl[ii];
	ii++;
	B243 = ictrl[ii];
	ii++;
    B260 = dctrl[dd];
	dd++;
    B264 = dctrl[dd];
	dd++;
    B267 = dctrl[dd];
	dd++;
    B269 = ictrl[ii];
	ii++;
	B270 = ictrl[ii];
	ii++;
    B274 = ictrl[ii];
	ii++;
    B281 = ictrl[ii];
	ii++;
    B291 = ictrl[ii];
	ii++;
    B295 = ictrl[ii];
	ii++;
    B308 = ictrl[ii];
	ii++;
    B309 = dctrl[dd];
	dd++;
    B310 = ictrl[ii];
	ii++;
	
	C1 = dctrl[dd];
	dd++;
    C2 = dctrl[dd];
	dd++;
	C3 = dctrl[dd];
	dd++;
	C4 = dctrl[dd];
	dd++;
	C5 = dctrl[dd];
	dd++;
    C10 = ictrl[ii];
	ii++;
	C15 = ictrl[ii];
	ii++;
	C20 = ictrl[ii];
	ii++;
    C50_1=dctrl[dd];
	dd++;
    C50_2=dctrl[dd];
	dd++;
    C51 = dctrl[dd];
	dd++;
    C52 = dctrl[dd];
	dd++;
    C53 = dctrl[dd];
	dd++;
    C54 = dctrl[dd];
	dd++;
    C55 = dctrl[dd];
	dd++;
    C56 = dctrl[dd];
	dd++;
    C57_1 = dctrl[dd];
	dd++;
    C57_2 = dctrl[dd];
	dd++;
    C57_3 = dctrl[dd];
	dd++;
    C57_4 = dctrl[dd];
	dd++;
    C58_1 = dctrl[dd];
	dd++;
    C58_2 = dctrl[dd];
	dd++;
    C58_3 = dctrl[dd];
	dd++;
    C58_4 = dctrl[dd];
	dd++;
	C75 = ictrl[ii];
	ii++;

    D10 = ictrl[ii];
	ii++;
	D11 = ictrl[ii];
	ii++;
    D20 = ictrl[ii];
	ii++;
	D21 = ictrl[ii];
	ii++;
	D22 = ictrl[ii];
	ii++;
	D23 = ictrl[ii];
	ii++;
    D24 = ictrl[ii];
	ii++;
	D29 = dctrl[dd];
	dd++;
    D30 = ictrl[ii];
	ii++;
    D31 = ictrl[ii];
	ii++;
	D32 = ictrl[ii];
	ii++;
	D33 = ictrl[ii];
	ii++;
	D34 = ictrl[ii];
	ii++;
	D35 = dctrl[dd];
	dd++;
    D36 = ictrl[ii];
	ii++;
    D38 = ictrl[ii];
	ii++;
	
	
    F10 = ictrl[ii];
	ii++;
    F11 = ictrl[ii];
	ii++;
	F19 = dctrl[dd];
	dd++;
    F30 = ictrl[ii];
	ii++;
    F31 = ictrl[ii];
	ii++;
    F32 = ictrl[ii];
	ii++;
    F33 = dctrl[dd];
	dd++;
    F34 = ictrl[ii];
	ii++;
    F35 = ictrl[ii];
	ii++;
	F36 = ictrl[ii];
	ii++;
	F39 = dctrl[dd];
	dd++;
    F40 = ictrl[ii];
	ii++;
    F41 = ictrl[ii];
	ii++;
	F42 = dctrl[dd];
	dd++;
    F43 = dctrl[dd];
	dd++;
    F44 = ictrl[ii];
	ii++;
    F45 = dctrl[dd];
	dd++;
    F46 = ictrl[ii];
	ii++;
    F47 = ictrl[ii];
	ii++;
	F49 = ictrl[ii];
	ii++;
    F50 = ictrl[ii];
	ii++;
    F50_flag = ictrl[ii];
	ii++;
    F51 = dctrl[dd];
	dd++;
    F52 = dctrl[dd];
	dd++;
    F53 = dctrl[dd];
	dd++;
    F54 = dctrl[dd];
	dd++;
    F55 = dctrl[dd];
	dd++;
    F56 = dctrl[dd];
	dd++;
    F57_1 = dctrl[dd];
	dd++;
    F57_2 = dctrl[dd];
	dd++;
    F57_3 = dctrl[dd];
	dd++;
    F57_4 = dctrl[dd];
	dd++;
    F58_1 = dctrl[dd];
	dd++;
    F58_2 = dctrl[dd];
	dd++;
    F58_3 = dctrl[dd];
	dd++;
    F58_4 = dctrl[dd];
	dd++;
    F59_xm = dctrl[dd];
	dd++;
    F59_ym = dctrl[dd];
	dd++;
    F59_zs = dctrl[dd];
	dd++;
    F59_ze = dctrl[dd];
	dd++;
    F59_r = dctrl[dd];
	dd++;
    F60 = dctrl[dd];
	dd++;
    F61 = dctrl[dd];
	dd++;
    F62 = dctrl[dd];
	dd++;
    F63 = dctrl[dd];
	dd++;
	F64 = ictrl[ii];
	ii++;
	F64_xs = dctrl[dd];
	dd++;
	F64_ys = dctrl[dd];
	dd++;
	F64_zs = dctrl[dd];
	dd++;
	F64_alpha = dctrl[dd];
	dd++;
    F70 = ictrl[ii];
	ii++;
	F71 = ictrl[ii];
	ii++;
	F72 = ictrl[ii];
	ii++;
    F80 = ictrl[ii];
	ii++;
    F84 = dctrl[dd];
	dd++;
    F85 = ictrl[ii];
	ii++;
    F101 = ictrl[ii];
	ii++;
	F150 = ictrl[ii];
	ii++;
	F151 = ictrl[ii];
	ii++;
	F300 = ictrl[ii];
	ii++;
	F305 = ictrl[ii];
	ii++;
	F310 = ictrl[ii];
	ii++;
	F321 = dctrl[dd];
	dd++;
	F322 = dctrl[dd];
	dd++;
	F323 = dctrl[dd];
	dd++;
	F350 = ictrl[ii];
	ii++;
	F360 = dctrl[dd];
	dd++;
	F361 = dctrl[dd];
	dd++;
	F362 = dctrl[dd];
	dd++;
    F369 = ictrl[ii];
	ii++;
	F370 = ictrl[ii];
	ii++;
	F371 = ictrl[ii];
	ii++;
    F374 = ictrl[ii];
	ii++;
    F375 = ictrl[ii];
	ii++;
    F378 = ictrl[ii];
	ii++;
    F379 = ictrl[ii];
	ii++;
	F380 = dctrl[dd];
	dd++;
	F381 = dctrl[dd];
	dd++;
	F382 = dctrl[dd];
	dd++;
	F390 = ictrl[ii];
	ii++;
	F391 = ictrl[ii];
	ii++;
    F394 = ictrl[ii];
	ii++;
    F395 = ictrl[ii];
	ii++;
    F398 = ictrl[ii];
	ii++;
    F399 = ictrl[ii];
	ii++;
	
	

    G1  = ictrl[ii];
	ii++;
    G10 = ictrl[ii];
	ii++;
    G11 = ictrl[ii];
	ii++;
    G12 = ictrl[ii];
	ii++;
    G20 = ictrl[ii];
	ii++;
    G21 = ictrl[ii];
	ii++;
    G22 = ictrl[ii];
	ii++;
    G30 = ictrl[ii];
	ii++;
    G40 = ictrl[ii];
	ii++;
    G50 = ictrl[ii];
	ii++;
    G52 = ictrl[ii];
	ii++;
    G60 = ictrl[ii];
	ii++;
    G61 = ictrl[ii];
	ii++;
	G81 = ictrl[ii];
	ii++;
	G95 = ictrl[ii];
	ii++;

    H1 = dctrl[dd];
	dd++;
    H2 = dctrl[dd];
	dd++;
    H10 = ictrl[ii];
	ii++;
    H50_1=dctrl[dd];
	dd++;
    H50_2=dctrl[dd];
	dd++;
    H51 = dctrl[dd];
	dd++;
    H52 = dctrl[dd];
	dd++;
    H53 = dctrl[dd];
	dd++;
    H54 = dctrl[dd];
	dd++;
    H55 = dctrl[dd];
	dd++;
    H56 = dctrl[dd];
	dd++;
    H57_1 = dctrl[dd];
	dd++;
    H57_2 = dctrl[dd];
	dd++;
    H57_3 = dctrl[dd];
	dd++;
    H57_4 = dctrl[dd];
	dd++;
    H58_1 = dctrl[dd];
	dd++;
    H58_2 = dctrl[dd];
	dd++;
    H58_3 = dctrl[dd];
	dd++;
    H58_4 = dctrl[dd];
	dd++;

    N5  = ictrl[ii];
    ii++;
    N8  = ictrl[ii];
    ii++;
    N9  = ictrl[ii];
    ii++;
	N10 = ictrl[ii];
    ii++;
    N11 = ictrl[ii];
    ii++;
    N12 = ictrl[ii];
    ii++;
    N13 = ictrl[ii];
    ii++;
    N14 = ictrl[ii];
    ii++;
	N15 = ictrl[ii];
    ii++;
	N16 = ictrl[ii];
    ii++;
	N17 = dctrl[dd];
    dd++;
	N18 = dctrl[dd];
    dd++;
    N40 = ictrl[ii];
    ii++;
    N41 = dctrl[dd];
    dd++;
	N42 = ictrl[ii];
    ii++;
	N43 = dctrl[dd];
    dd++;
    N44 = dctrl[dd];
    dd++;
    N45 = ictrl[ii];
    ii++;
    N46 = ictrl[ii];
    ii++;
    N47 = dctrl[dd];
    dd++;
    N48 = ictrl[ii];
    ii++;
    N49 = dctrl[dd];
    dd++;
    N50 = ictrl[ii];
    ii++;
    N51 = dctrl[dd];
    dd++;
    N52 = ictrl[ii];
    ii++;
    N53 = dctrl[dd];
    dd++;
    N54 = dctrl[dd];
    dd++;
    N55 = dctrl[dd];
    dd++;
    N56 = dctrl[dd];
    dd++;
    N57_1 = ictrl[ii];
    ii++;
    N57_2 = ictrl[ii];
    ii++;
    N57_3 = dctrl[dd];
    dd++;
    N58 = ictrl[ii];
    ii++;
    N60 = ictrl[ii];
    ii++;
    N61 = dctrl[dd];
    dd++;
    N71 = ictrl[ii];
    ii++;
    N72 = ictrl[ii];
    ii++;
    N73 = ictrl[ii];
    ii++;
	

    P10 = ictrl[ii];
	ii++;
	P11 = ictrl[ii];
	ii++;
    P12 = ictrl[ii];
	ii++;
	P14 = ictrl[ii];
	ii++;
    P15 = ictrl[ii];
	ii++;
    P17 = ictrl[ii];
	ii++;
	P18 = ictrl[ii];
	ii++;
	P19 = ictrl[ii];
	ii++;
	P20 = ictrl[ii];
	ii++;
	P23 = ictrl[ii];
	ii++;
    P24 = ictrl[ii];
	ii++;
    P25 = ictrl[ii];
	ii++;
	P26 = ictrl[ii];
	ii++;
	P27 = ictrl[ii];
	ii++;
	P28 = ictrl[ii];
	ii++;
	P29 = ictrl[ii];
	ii++;
    P30 = dctrl[dd];
	dd++;
	P34 = dctrl[dd];
	dd++;
    P35 = ictrl[ii];
	ii++;
    P40 = ictrl[ii];
	ii++;
    P41 = ictrl[ii];
	ii++;
	P42 = dctrl[dd];
	dd++;
    P50 = ictrl[ii];
	ii++;
	P51 = ictrl[ii];
	ii++;
    P52 = ictrl[ii];
	ii++;
    P53 = ictrl[ii];
	ii++;
	P54 = ictrl[ii];
	ii++;
	P55 = dctrl[dd];
	dd++;
	P56 = ictrl[ii];
	ii++;
	P59 = ictrl[ii];
	ii++;
	P61 = ictrl[ii];
	ii++;
	P62 = ictrl[ii];
	ii++;
	P66 = ictrl[ii];
	ii++;
	P67 = ictrl[ii];
	ii++;
    P71 = ictrl[ii];
	ii++;
    P75 = ictrl[ii];
	ii++;
	P78 = ictrl[ii];
	ii++;
    P79 = ictrl[ii];
	ii++;
    P81 = ictrl[ii];
	ii++;
    P82_x = dctrl[dd];
	dd++;
    P82_y = dctrl[dd];
	dd++;
	P83 = dctrl[dd];
	dd++;
	P84 = dctrl[dd];
	dd++;
	P85 = ictrl[ii];
	ii++;
    P86_x = dctrl[dd];
	dd++;
    P86_y = dctrl[dd];
	dd++;
	P87 = dctrl[dd];
	dd++;
	P88 = dctrl[dd];
	dd++;
	P89_cm = dctrl[dd];
	dd++;
	P89_cd = dctrl[dd];
	dd++;
	P90 = dctrl[dd];
	dd++;
	P91 = dctrl[dd];
	dd++;
    P92 = ictrl[ii];
	ii++;
    P93 = ictrl[ii];
	ii++;
	P101 = ictrl[ii];
	ii++;
	P101_xm = dctrl[dd];
	dd++;
	P101_ym = dctrl[dd];
	dd++;
	P101_zs = dctrl[dd];
	dd++;
	P101_ze = dctrl[dd];
	dd++;
	P101_r1 = dctrl[dd];
	dd++;
	P101_r2 = dctrl[dd];
	dd++;
	P121 = ictrl[ii];
	ii++;
	P122 = ictrl[ii];
	ii++;
	P123 = ictrl[ii];
	ii++;
	P124 = ictrl[ii];
	ii++;
	P125 = ictrl[ii];
	ii++;
	P126 = ictrl[ii];
	ii++;
	P151 = ictrl[ii];
	ii++;
	P152 = ictrl[ii];
	ii++;
	P180 = ictrl[ii];
	ii++;
	P181 = ictrl[ii];
	ii++;
	P182 = dctrl[dd];
	dd++;
    P210 = ictrl[ii];
	ii++;
	P211 = ictrl[ii];
	ii++;
	P212 = dctrl[dd];
	dd++;
    P230 = ictrl[ii];
	ii++;
    P240 = ictrl[ii];
	ii++;
	P351 = ictrl[ii];
	ii++;
	P352 = ictrl[ii];
	ii++;


    S10 = ictrl[ii];
	ii++;
    S11 = ictrl[ii];
	ii++;
    S12 = ictrl[ii];
	ii++;
    S13 = dctrl[dd];
	dd++;
    S14 = dctrl[dd];
	dd++;
    S15 = ictrl[ii];
	ii++;
    S16 = ictrl[ii];
	ii++;
    S17 = ictrl[ii];
	ii++;
    S18 = ictrl[ii];
	ii++;
	S19 = dctrl[dd];
	dd++;
    S20 = dctrl[dd];
	dd++;
	S21 = dctrl[dd];
	dd++;
    S22 = dctrl[dd];
	dd++;
    S23 = dctrl[dd];
	dd++;
    S24 = dctrl[dd];
	dd++;
    S25 = dctrl[dd];
	dd++;
	S28 = dctrl[dd];
	dd++;
	S29 = dctrl[dd];
	dd++;
    S30 = dctrl[dd];
	dd++;
	S31 = ictrl[ii];
	ii++;
    S37 = ictrl[ii];
	ii++;
	S38 = ictrl[ii];
	ii++;
	S39 = ictrl[ii];
	ii++;
    S41 = ictrl[ii];
	ii++;
	S42 = ictrl[ii];
	ii++;
    S43 = ictrl[ii];
	ii++;
    S44 = ictrl[ii];
	ii++;
	S45 = dctrl[dd];
	dd++;
	S46 = dctrl[dd];
	dd++;
	S47 = dctrl[dd];
	dd++;
	S48 = dctrl[dd];
	dd++;
    S50 = ictrl[ii];
	ii++;
    S57 = dctrl[dd];
	dd++;
    S60 = dctrl[dd];
	dd++;
    S71 = dctrl[dd];
	dd++;
    S72 = dctrl[dd];
	dd++;
    S73 = ictrl[ii];
	ii++;
    S77 = ictrl[ii];
	ii++;
    S77_xs = dctrl[dd];
	dd++;
    S77_xe = dctrl[dd];
	dd++;
    S80 = ictrl[ii];
	ii++;
    S81 = dctrl[dd];
	dd++;
    S82 = dctrl[dd];
	dd++;
    S90 = ictrl[ii];
	ii++;
    S91 = ictrl[ii];
	ii++;
	S92 = dctrl[dd];
	dd++;
    S93 = dctrl[dd];
	dd++;
	S100 = ictrl[ii];
	ii++;
	S101 = ictrl[ii];
	ii++;
    S102 = ictrl[ii];
	ii++;
    S103 = ictrl[ii];
	ii++;
    S116 = dctrl[dd];
	dd++;
    S117 = dctrl[dd];
	dd++;

    T10 = ictrl[ii];
    ii++;
    T11 = ictrl[ii];
    ii++;
    T12 = ictrl[ii];
    ii++;
	T13 = dctrl[dd];
    dd++;
    T30 = ictrl[ii];
    ii++;
    T31 = dctrl[dd];
    dd++;
	T35 = dctrl[dd];
    dd++;
	T36 = ictrl[ii];
    ii++;
	T37 = dctrl[dd];
    dd++;
    T38 = dctrl[dd];
    dd++;
    T39 = dctrl[dd];
    dd++;
    T40 = ictrl[ii];
    ii++;
	T41 = ictrl[ii];
    ii++;
	T42 = dctrl[dd];
    dd++;
	T43 = ictrl[ii];
    ii++;
    T51 = dctrl[dd];
    dd++;
	T52 = dctrl[dd];
    dd++;
	T53 = dctrl[dd];
    dd++;
	
    W1  = dctrl[dd];
    dd++;
    W2  = dctrl[dd];
    dd++;
    W3  = dctrl[dd];
    dd++;
    W4  = dctrl[dd];
    dd++;
    W5  = dctrl[dd];
    dd++;
	W6  = dctrl[dd];
    dd++;
	W7  = dctrl[dd];
    dd++;
    W10 = dctrl[dd];
    dd++;
    W20 = dctrl[dd];
    dd++;
    W21 = dctrl[dd];
    dd++;
    W22 = dctrl[dd];
    dd++;
	W30 = ictrl[ii];
    ii++;
	W31 = dctrl[dd];
    dd++;
    W90 = ictrl[ii];
    ii++;
	W95 = dctrl[dd];
    dd++;
    W96 = dctrl[dd];
    dd++;
    W97 = dctrl[dd];
    dd++;
	W98 = dctrl[dd];
    dd++;
    W101 = ictrl[ii];
    ii++;
    W102_phi = dctrl[dd];
    dd++;
    W102_c = dctrl[dd];
    dd++;
    W103 = dctrl[dd];
    dd++;
    W104 = dctrl[dd];
    dd++;
    W110 = ictrl[ii];
    ii++;
    W111 = ictrl[ii];
    ii++;
    W112 = dctrl[dd];
    dd++;
	
	X10 = ictrl[ii];
	ii++;
	X11_u = ictrl[ii];
	ii++;
	X11_v = ictrl[ii];
	ii++;
	X11_w = ictrl[ii];
	ii++;
	X11_p = ictrl[ii];
	ii++;
	X11_q = ictrl[ii];
	ii++;
	X11_r = ictrl[ii];
	ii++;
    X12 = ictrl[ii];
	ii++;
	X13 = ictrl[ii];
	ii++;
	X18 = ictrl[ii];
	ii++;
	X19 = ictrl[ii];
	ii++;
	X21 = ictrl[ii];
	ii++;
	X21_d = dctrl[dd];
	dd++;
	X22 = ictrl[ii];
	ii++;
	X22_m = dctrl[dd];
	dd++;
	X23 = ictrl[ii];
	ii++;
	X23_x = dctrl[dd];
	dd++;
	X23_y = dctrl[dd];
	dd++;
	X23_z = dctrl[dd];
	dd++;
	X24 = ictrl[ii];
	ii++;
	X24_Ix = dctrl[dd];
	dd++;
	X24_Iy = dctrl[dd];
	dd++;
	X24_Iz = dctrl[dd];
	dd++;
	X25_Cp = dctrl[dd];
	dd++;
	X25_Cq = dctrl[dd];
	dd++;
	X25_Cr = dctrl[dd];
	dd++;
    X26_Ku = dctrl[dd];
	dd++;
	X26_Kv = dctrl[dd];
	dd++;
	X26_Kw = dctrl[dd];
	dd++;
	X27 = ictrl[ii];
	ii++;
    X27_x = dctrl[dd];
	dd++;
    X27_y = dctrl[dd];
	dd++;
    X27_z = dctrl[dd];
	dd++;
    X31 = ictrl[ii];
	ii++;
	X32 = ictrl[ii];
	ii++;
	X33 = ictrl[ii];
	ii++;
    X34 = ictrl[ii];
	ii++;
	X38 = ictrl[ii];
	ii++;
    X40 = ictrl[ii];
	ii++;
	X41 = dctrl[dd];
	dd++;
	X100 = ictrl[ii];
	ii++;
	X100_x = dctrl[dd];
	dd++;
	X100_y = dctrl[dd];
	dd++;
	X100_z = dctrl[dd];
	dd++;
	X101 = ictrl[ii];
	ii++;
	X101_phi = dctrl[dd];
	dd++;
	X101_theta = dctrl[dd];
	dd++;
	X101_psi = dctrl[dd];
	dd++;
	X102 = ictrl[ii];
	ii++;
	X102_u = dctrl[dd];
	dd++;
	X102_v = dctrl[dd];
	dd++;
	X102_w = dctrl[dd];
	dd++;
	X103 = ictrl[ii];
	ii++;
	X103_p = dctrl[dd];
	dd++;
	X103_q = dctrl[dd];
	dd++;
	X103_r = dctrl[dd];
	dd++;
	X110 = ictrl[ii];
	ii++;
	X120 = ictrl[ii];
	ii++;
	X120_rad = dctrl[dd];
	dd++;
	X120_xc = dctrl[dd];
	dd++;
	X120_yc = dctrl[dd];
	dd++;
	X120_zc = dctrl[dd];
	dd++;
	X131 = ictrl[ii];
	ii++;
	X131_rad = dctrl[dd];
	dd++;
	X131_h = dctrl[dd];
	dd++;
	X131_xc = dctrl[dd];
	dd++;
	X131_yc = dctrl[dd];
	dd++;
	X131_zc = dctrl[dd];
	dd++;
	X132 = ictrl[ii];
	ii++;
	X132_rad = dctrl[dd];
	dd++;
	X132_h = dctrl[dd];
	dd++;
	X132_xc = dctrl[dd];
	dd++;
	X132_yc = dctrl[dd];
	dd++;
	X132_zc = dctrl[dd];
	dd++;
	X133 = ictrl[ii];
	ii++;
	X133_rad = dctrl[dd];
	dd++;
	X133_h = dctrl[dd];
	dd++;
	X133_xc = dctrl[dd];
	dd++;
	X133_yc = dctrl[dd];
	dd++;
	X133_zc = dctrl[dd];
	dd++;
	X153 = ictrl[ii];
	ii++;
	X153_xs = dctrl[dd];
	dd++;
	X153_xe = dctrl[dd];
	dd++;
	X153_ys = dctrl[dd];
	dd++;
	X153_ye = dctrl[dd];
	dd++;
	X153_zs = dctrl[dd];
	dd++;
	X153_ze = dctrl[dd];
	dd++;
    X163 = ictrl[ii];
	ii++;
    X164 = ictrl[ii];
	ii++;
    X180 = ictrl[ii];
	ii++;
    X181 = dctrl[dd];
	dd++;
    X182 = ictrl[ii];
	ii++;
    X182_x = dctrl[dd];
	dd++;
    X182_y = dctrl[dd];
	dd++;
    X182_z = dctrl[dd];
	dd++;
    X183 = ictrl[ii];
	ii++;
    X183_x = dctrl[dd];
	dd++;
    X183_y = dctrl[dd];
	dd++;
    X183_z = dctrl[dd];
	dd++;
    X183_phi = dctrl[dd];
	dd++;
    X183_theta = dctrl[dd];
	dd++;
    X183_psi = dctrl[dd];
	dd++;
	X210 = ictrl[ii];
	ii++;
	X210_u = dctrl[dd];
	dd++;
	X210_v = dctrl[dd];
	dd++;
	X210_w = dctrl[dd];
	dd++;
	X211 = ictrl[ii];
	ii++;
	X211_p = dctrl[dd];
	dd++;
	X211_q = dctrl[dd];
	dd++;
	X211_r = dctrl[dd];
	dd++;	
    X221 = ictrl[ii];
	ii++;
	X221_xs = dctrl[dd];
	dd++;
    X221_xe = dctrl[dd];
	dd++;
    X221_ys = dctrl[dd];
	dd++;
    X221_ye = dctrl[dd];
	dd++;
    X221_zs = dctrl[dd];
	dd++;
    X221_ze = dctrl[dd];
	dd++;
	X310 = ictrl[ii];
	ii++;	
    X311 = ictrl[ii];
	ii++;
    X320 = ictrl[ii];
	ii++;
	
	Y50 = ictrl[ii];
	ii++;
	Y60 = ictrl[ii];
	ii++;
    Y71 = ictrl[ii];
	ii++;
    Y72 = ictrl[ii];
	ii++;
    Y73 = ictrl[ii];
	ii++;
    Y74 = ictrl[ii];
	ii++;
    
    

// --------------------------	
	
	if(B67>0)
	{
	Darray(B67_val,B67);
	Darray(B67_dist,B67);
	Darray(B67_b,B67);
	Darray(B67_x,B67);
	Darray(B67_y,B67);
	}
	
	if(B68>0)
	{
	Darray(B68_val,B68);
	Darray(B68_dist,B68);
	Darray(B68_b,B68);
	Darray(B68_x,B68);
	Darray(B68_y,B68);
	}
	
	if(B69>0)
	{
	Darray(B69_val,B69);
	Darray(B69_dist,B69);
	Darray(B69_b,B69);
	Darray(B69_x,B69);
	Darray(B69_y,B69);
	}
	
	if(B70>0)
	{
	Darray(B70_val,B70);
	Darray(B70_dist,B70);
	Darray(B70_b,B70);
	Darray(B70_x,B70);
	Darray(B70_y,B70);
	}
	
	if(B71>0)
	{
	Darray(B71_val,B71);
	Darray(B71_dist,B71);
	Darray(B71_b,B71);
	Darray(B71_x,B71);
	Darray(B71_y,B71);
	}
	
	if(B106>0)
	{
	Darray(B106_b,B106);
	Darray(B106_x,B106);
	Darray(B106_y,B106);
	}
	
	if(B107>0)
	{
	Darray(B107_xs,B107);
	Darray(B107_xe,B107);
	Darray(B107_ys,B107);
    Darray(B107_ye,B107);
    Darray(B107_d,B107);
	}
    
    if(B108>0)
	{
	Darray(B108_xs,B108);
	Darray(B108_xe,B108);
	Darray(B108_ys,B108);
    Darray(B108_ye,B108);
    Darray(B108_d,B108);
	}
    
	if(B240>0)
	{	
	Darray(B240_C,B240);
    Darray(B240_D,B240);
    Darray(B240_xs,B240);
	Darray(B240_xe,B240);
	Darray(B240_ys,B240);
	Darray(B240_ye,B240);
	Darray(B240_zs,B240);
	Darray(B240_ze,B240);
	}
    
	if(B270>0)
	{	
	Darray(B270_xs,B270);
	Darray(B270_xe,B270);
	Darray(B270_ys,B270);
	Darray(B270_ye,B270);
	Darray(B270_zs,B270);
	Darray(B270_ze,B270);
    Darray(B270_n,B270);
    Darray(B270_d50,B270);
	Darray(B270_alpha,B270);
	Darray(B270_beta,B270);
	}
    
    if(B274>0)
	{	
	Darray(B274_xc,B274);
	Darray(B274_yc,B274);
	Darray(B274_zs,B274);
	Darray(B274_ze,B274);
	Darray(B274_r,B274);
    Darray(B274_n,B274);
    Darray(B274_d50,B274);
	Darray(B274_alpha,B274);
	Darray(B274_beta,B274);
	}
    
    if(B281>0)
	{	
	Darray(B281_xs,B281);
	Darray(B281_xe,B281);
	Darray(B281_ys,B281);
	Darray(B281_ye,B281);
	Darray(B281_zs,B281);
	Darray(B281_ze,B281);
    Darray(B281_n,B281);
    Darray(B281_d50,B281);
	Darray(B281_alpha,B281);
	Darray(B281_beta,B281);
	}
    
    if(B291>0)
	{	
	Darray(B291_xs,B291);
	Darray(B291_xe,B291);
	Darray(B291_ys,B291);
	Darray(B291_ye,B291);
	Darray(B291_zs,B291);
	Darray(B291_ze,B291);
    Darray(B291_d,B291);
    Darray(B291_n,B291);
    Darray(B291_d50,B291);
	Darray(B291_alpha,B291);
	Darray(B291_beta,B291);
	}
    
    if(B310>0)
	{	
	Darray(B310_xs,B310);
	Darray(B310_xe,B310);
	Darray(B310_ys,B310);
	Darray(B310_ye,B310);
	Darray(B310_zs,B310);
	Darray(B310_ze,B310);
    Darray(B310_N,B310);
    Darray(B310_D,B310);
	Darray(B310_Cd,B310);
	}
	
	if(C75>0)
	{
	Darray(C75_x,C75);   
	Darray(C75_z,C75);  
	
	Darray(C75_a,C75);  	
	Darray(C75_s,C75);  
	Darray(C75_l,C75);  
	Darray(C75_v,C75);  
	}

	if(F70>0)
	{
	Darray(F70_xs,F70);  
	Darray(F70_xe,F70);  
	
	Darray(F70_ys,F70);  
	Darray(F70_ye,F70);  
	
	Darray(F70_zs,F70);  
	Darray(F70_ze,F70);  
	}
	
	if(F71>0)
	{
	Darray(F71_xs,F71);  
	Darray(F71_xe,F71);  
	
	Darray(F71_ys,F71);  
	Darray(F71_ye,F71);  
	
	Darray(F71_zs,F71);  
	Darray(F71_ze,F71);  
	}
	
	if(F72>0)
	{
	Darray(F72_xs,F72);  
	Darray(F72_xe,F72);  
	
	Darray(F72_ys,F72);  
	Darray(F72_ye,F72);  
	
	Darray(F72_h,F72);  
	}
	if(F369>0)
	{
	Darray(F369_x,F369);   
	Darray(F369_z,F369);  
	
	Darray(F369_a,F369);  	
	Darray(F369_s,F369);  
	Darray(F369_l,F369);  
	Darray(F369_v,F369);  
	}
    
	if(F370>0)
	{
	Darray(F370_xs,F370);  
	Darray(F370_xe,F370);  
	
	Darray(F370_ys,F370);  
	Darray(F370_ye,F370);  
	
	Darray(F370_zs,F370);  
	Darray(F370_ze,F370);  
	}
	
	if(F371>0)
	{
	Darray(F371_xs,F371);  
	Darray(F371_xe,F371);  
	
	Darray(F371_ys,F371);  
	Darray(F371_ye,F371);  
	
	Darray(F371_zs,F371);  
	Darray(F371_ze,F371);  
	}
	
    if(F374>0)
	{
	Darray(F374_xc,F374);  
	Darray(F374_zc,F374);
    Darray(F374_r,F374);  
	}
    
    if(F375>0)
	{
	Darray(F375_xc,F375);  
	Darray(F375_zc,F375);
    Darray(F375_r,F375);  
	}
    
    if(F378>0)
	{
	Darray(F378_xc,F378);  
    Darray(F378_yc,F378);  
	Darray(F378_zc,F378);
    Darray(F378_r,F378);  
	}
	
    if(F379>0)
	{
	Darray(F379_xc,F379);  
    Darray(F379_yc,F379);  
	Darray(F379_zc,F379);
    Darray(F379_r,F379);  
	}
	
	if(F390>0)
	{
	Darray(F390_xs,F390);  
	Darray(F390_xe,F390);  
	
	Darray(F390_ys,F390);  
	Darray(F390_ye,F390);  
	
	Darray(F390_zs,F390);  
	Darray(F390_ze,F390);  
	}
	
	if(F391>0)
	{
	Darray(F391_xs,F391);  
	Darray(F391_xe,F391);  
	
	Darray(F391_ys,F391);  
	Darray(F391_ye,F391);  
	
	Darray(F391_zs,F391);  
	Darray(F391_ze,F391);  
	}
    
    if(F394>0)
	{
	Darray(F394_xc,F394);  
	Darray(F394_zc,F394);
    Darray(F394_r,F394);  
	}
    
    if(F395>0)
	{
	Darray(F395_xc,F395);  
	Darray(F395_zc,F395);
    Darray(F395_r,F395);  
	}
    
    if(F398>0)
	{
	Darray(F398_xc,F398);  
    Darray(F398_yc,F398);  
	Darray(F398_zc,F398);
    Darray(F398_r,F398);  
	}
	
    if(F399>0)
	{
	Darray(F399_xc,F399);  
    Darray(F399_yc,F399);  
	Darray(F399_zc,F399);
    Darray(F399_r,F399);  
	}
	
	if(G60>0)
	{
	Darray(G60_xs,G60); 
	Darray(G60_xe,G60); 
	
	Darray(G60_ys,G60); 
	Darray(G60_ye,G60); 
	
	Darray(G60_zs,G60); 
	Darray(G60_ze,G60); 
	}
	
	if(G61>0)
	{
	Darray(G61_xs,G61); 
	Darray(G61_xe,G61); 
	
	Darray(G61_ys,G61); 
	Darray(G61_ye,G61); 
	
	Darray(G61_zs,G61); 
	Darray(G61_ze,G61); 
	}
	
	if(G81>0)
	{
	Darray(G81_xs,G81); 
	Darray(G81_xe,G81); 
	
	Darray(G81_ys,G81); 
	Darray(G81_ye,G81); 
	
	Darray(G81_zs,G81); 
	Darray(G81_ze,G81); 
	}
	
	if(G95>0)
	{
	Darray(G95_xs,G95); 
	Darray(G95_xe,G95); 
	
	Darray(G95_ys,G95); 
	Darray(G95_ye,G95); 
	
	Darray(G95_zs,G95); 
	Darray(G95_ze,G95); 
	}

	if(P35>0)
	{
    Darray(P35_ts,P35);  
	Darray(P35_te,P35);  
	Darray(P35_dt,P35);  
	}
	
	if(P50>0)
	{
    Darray(P50_x,P50);  
	Darray(P50_y,P50);  
	}
	
	if(P51>0)
	{
    Darray(P51_x,P51);  
	Darray(P51_y,P51);  
	}
	
	if(P52>0)
	Darray(P52_y,P52);

	if(P56>0)
	Darray(P56_x,P56);  
	
	if(P61>0)
	{
    Darray(P61_x,P61);  
	Darray(P61_y,P61); 
	Darray(P61_z,P61);  
	}
	
	if(P62>0)
	{
	Darray(P62_xs,P62); 
	Darray(P62_xe,P62); 
	
	Darray(P62_ys,P62); 
	Darray(P62_ye,P62); 
	
	Darray(P62_zs,P62); 
	Darray(P62_ze,P62); 
	}
	
	if(P67>0)
	Darray(P67_x,P67);  
	
	if(P81>0)
	{
	Darray(P81_xs,P81); 
	Darray(P81_xe,P81); 
	
	Darray(P81_ys,P81); 
	Darray(P81_ye,P81); 
	
	Darray(P81_zs,P81); 
	Darray(P81_ze,P81); 
	}
	
	if(P85>0)
	{
	Darray(P85_xs,P85); 
	Darray(P85_xe,P85); 
	
	Darray(P85_ys,P85); 
	Darray(P85_ye,P85); 
	
	Darray(P85_zs,P85); 
	Darray(P85_ze,P85); 
	}
	
	if(P121>0)
	{
    Darray(P121_x,P121);  
	Darray(P121_y,P121);  
	}
	
	if(P123>0)
	Darray(P123_y,P123);
	
	if(P124>0)
	Darray(P124_x,P124);
	
	if(P125>0)
	{
    Darray(P125_x,P125);  
	Darray(P125_y,P125);  
	}
    
    if(P230>0)
	{
    Darray(P230_x,P230);  
	}
    
    if(P240>0)
	{
    Darray(P240_x,P240);  
	}
	
	if(P351>0)
	{
    Darray(P351_x,P351);  
	Darray(P351_y,P351);  
	}
	
	if(P352>0)
	{
    Darray(P352_x,P352);  
	Darray(P352_y,P352);  
	}
    
    if(S73>0)
	{
	Darray(S73_val,S73);
	Darray(S73_dist,S73);
	Darray(S73_b,S73);
	Darray(S73_x,S73);
	Darray(S73_y,S73);
	}
	
	if(X110>0)
	{
	Darray(X110_xs,X110);  
	Darray(X110_xe,X110);  
	
	Darray(X110_ys,X110);  
	Darray(X110_ye,X110);  
	
	Darray(X110_zs,X110);  
	Darray(X110_ze,X110);  
	}
    
    if(X163>0)
	{
	Darray(X163_x1,X163);
    Darray(X163_y1,X163);
    Darray(X163_z1,X163);
    Darray(X163_x2,X163);
    Darray(X163_y2,X163);
    Darray(X163_z2,X163);
    Darray(X163_x3,X163);
    Darray(X163_y3,X163);
    Darray(X163_z3,X163);
    Darray(X163_x4,X163);
    Darray(X163_y4,X163);
    Darray(X163_z4,X163);
    Darray(X163_x5,X163);
    Darray(X163_y5,X163);
    Darray(X163_z5,X163);
    Darray(X163_x6,X163);
    Darray(X163_y6,X163);
    Darray(X163_z6,X163);
    }
    
    if(X164>0)
	{
    Darray(X164_x1,X164);
    Darray(X164_y1,X164);
    Darray(X164_z1,X164);
    Darray(X164_x2,X164);
    Darray(X164_y2,X164);
    Darray(X164_z2,X164);
    Darray(X164_x3,X164);
    Darray(X164_y3,X164);
    Darray(X164_z3,X164);
    Darray(X164_x4,X164);
    Darray(X164_y4,X164);
    Darray(X164_z4,X164);
    Darray(X164_x5,X164);
    Darray(X164_y5,X164);
    Darray(X164_z5,X164);
    Darray(X164_x6,X164);
    Darray(X164_y6,X164);
    Darray(X164_z6,X164);
    Darray(X164_x7,X164);
    Darray(X164_y7,X164);
    Darray(X164_z7,X164);
    Darray(X164_x8,X164);
    Darray(X164_y8,X164);
    Darray(X164_z8,X164);
	}
    
    if(X311>0)
	{
		Darray(X311_xs,X311);  
		Darray(X311_xe,X311);  
		Darray(X311_ys,X311);  
		Darray(X311_ye,X311);  
		Darray(X311_zs,X311);  
		Darray(X311_ze,X311);  
	
		Darray(X311_w,X311); 
		Darray(X311_rho_c,X311); 
		Darray(X311_EA,X311); 
		Darray(X311_d,X311);  
		Darray(X311_l,X311); 
		Darray(X311_H,X311); 
		Darray(X311_P,X311); 
		Darray(X311_facT,X311);  
	}


// --------------------------

	
	for(n=0;n<B67;++n)
    {
    B67_val[n]= dctrl[dd];
    dd++;
	B67_dist[n]= dctrl[dd];
    dd++;
	B67_b[n]  = dctrl[dd];
    dd++;
    B67_x[n]  = dctrl[dd];
    dd++;
	B67_y[n]  = dctrl[dd];
    dd++;
    }
	
	for(n=0;n<B68;++n)
    {
    B68_val[n]= dctrl[dd];
    dd++;
	B68_dist[n]= dctrl[dd];
    dd++;
	B68_b[n]  = dctrl[dd];
    dd++;
    B68_x[n]  = dctrl[dd];
    dd++;
	B68_y[n]  = dctrl[dd];
    dd++;
    }
	
	for(n=0;n<B69;++n)
    {
    B69_val[n]= dctrl[dd];
    dd++;
	B69_dist[n]= dctrl[dd];
    dd++;
	B69_b[n]  = dctrl[dd];
    dd++;
    B69_x[n]  = dctrl[dd];
    dd++;
	B69_y[n]  = dctrl[dd];
    dd++;
    }
	
	for(n=0;n<B70;++n)
    {
    B70_val[n]= dctrl[dd];
    dd++;
	B70_dist[n]= dctrl[dd];
    dd++;
	B70_b[n]  = dctrl[dd];
    dd++;
    B70_x[n]  = dctrl[dd];
    dd++;
	B70_y[n]  = dctrl[dd];
    dd++;
    }
	
	for(n=0;n<B71;++n)
    {
    B71_val[n]= dctrl[dd];
    dd++;
	B71_dist[n]= dctrl[dd];
    dd++;
	B71_b[n]  = dctrl[dd];
    dd++;
    B71_x[n]  = dctrl[dd];
    dd++;
	B71_y[n]  = dctrl[dd];
    dd++;
    }
		
	for(n=0;n<B106;++n)
    {
    B106_b[n]  = dctrl[dd];
    dd++;
    B106_x[n]  = dctrl[dd];
    dd++;
	B106_y[n]  = dctrl[dd];
    dd++;
    }
	
	for(n=0;n<B107;++n)
    {
    B107_xs[n]  = dctrl[dd];
    dd++;
    B107_xe[n]  = dctrl[dd];
    dd++;
	B107_ys[n]  = dctrl[dd];
    dd++;
    B107_ye[n]  = dctrl[dd];
    dd++;
    B107_d[n]  = dctrl[dd];
    dd++;
    }
    
    for(n=0;n<B108;++n)
    {
    B108_xs[n]  = dctrl[dd];
    dd++;
    B108_xe[n]  = dctrl[dd];
    dd++;
	B108_ys[n]  = dctrl[dd];
    dd++;
    B108_ye[n]  = dctrl[dd];
    dd++;
    B108_d[n]  = dctrl[dd];
    dd++;
    }
	
	for(n=0;n<B240;++n)
    {
    B240_C[n]  = dctrl[dd];
    dd++;
    B240_D[n]  = dctrl[dd];
    dd++;
	B240_xs[n] = dctrl[dd];
    dd++;
    B240_xe[n] = dctrl[dd];
    dd++;
    B240_ys[n] = dctrl[dd];
    dd++;
    B240_ye[n] = dctrl[dd];
    dd++;
    B240_zs[n] = dctrl[dd];
    dd++;
    B240_ze[n] = dctrl[dd];
    dd++;
    }
    
    for(n=0;n<B270;++n)
    {
	B270_xs[n] = dctrl[dd];
    dd++;
    B270_xe[n] = dctrl[dd];
    dd++;
    B270_ys[n] = dctrl[dd];
    dd++;
    B270_ye[n] = dctrl[dd];
    dd++;
    B270_zs[n] = dctrl[dd];
    dd++;
    B270_ze[n] = dctrl[dd];
    dd++;
    B270_n[n]  = dctrl[dd];
    dd++;
    B270_d50[n]= dctrl[dd];
    dd++;
	B270_alpha[n]= dctrl[dd];
    dd++;
	B270_beta[n]= dctrl[dd];
    dd++;
    }
    
    for(n=0;n<B274;++n)
    {
	B274_xc[n] = dctrl[dd];
    dd++;
    B274_yc[n] = dctrl[dd];
    dd++;
    B274_zs[n] = dctrl[dd];
    dd++;
    B274_ze[n] = dctrl[dd];
    dd++;
    B274_r[n] = dctrl[dd];
    dd++;
    B274_n[n]  = dctrl[dd];
    dd++;
    B274_d50[n]= dctrl[dd];
    dd++;
	B274_alpha[n]= dctrl[dd];
    dd++;
	B274_beta[n]= dctrl[dd];
    dd++;
    }
    
    for(n=0;n<B281;++n)
    {
	B281_xs[n] = dctrl[dd];
    dd++;
    B281_xe[n] = dctrl[dd];
    dd++;
    B281_ys[n] = dctrl[dd];
    dd++;
    B281_ye[n] = dctrl[dd];
    dd++;
    B281_zs[n] = dctrl[dd];
    dd++;
    B281_ze[n] = dctrl[dd];
    dd++;
    B281_n[n]  = dctrl[dd];
    dd++;
    B281_d50[n]= dctrl[dd];
    dd++;
	B281_alpha[n]= dctrl[dd];
    dd++;
	B281_beta[n]= dctrl[dd];
    dd++;
    }
    
    for(n=0;n<B310;++n)
    {
	B310_xs[n] = dctrl[dd];
    dd++;
    B310_xe[n] = dctrl[dd];
    dd++;
    B310_ys[n] = dctrl[dd];
    dd++;
    B310_ye[n] = dctrl[dd];
    dd++;
    B310_zs[n] = dctrl[dd];
    dd++;
    B310_ze[n] = dctrl[dd];
    dd++;
    B310_N[n]  = dctrl[dd];
    dd++;
    B310_D[n]  = dctrl[dd];
    dd++;
	B310_Cd[n] = dctrl[dd];
    dd++;
    }
    
    for(n=0;n<B291;++n)
    {
	B291_xs[n] = dctrl[dd];
    dd++;
    B291_xe[n] = dctrl[dd];
    dd++;
    B291_ys[n] = dctrl[dd];
    dd++;
    B291_ye[n] = dctrl[dd];
    dd++;
    B291_zs[n] = dctrl[dd];
    dd++;
    B291_ze[n] = dctrl[dd];
    dd++;
    B291_d[n] = dctrl[dd];
    dd++;
    B291_n[n]  = dctrl[dd];
    dd++;
    B291_d50[n]= dctrl[dd];
    dd++;
	B291_alpha[n]= dctrl[dd];
    dd++;
	B291_beta[n]= dctrl[dd];
    dd++;
    }
	
	for(n=0;n<C75;++n)
    {
    C75_x[n] = dctrl[dd];
    dd++;
    C75_z[n] = dctrl[dd];
    dd++;
    C75_a[n] = dctrl[dd];
    dd++;
    C75_s[n] = dctrl[dd];
    dd++;
    C75_l[n] = dctrl[dd];
    dd++;
	C75_v[n] = dctrl[dd];
    dd++;
    }
	
    for(n=0;n<F70;++n)
    {
    F70_xs[n] = dctrl[dd];
    dd++;
    F70_xe[n] = dctrl[dd];
    dd++;
    F70_ys[n] = dctrl[dd];
    dd++;
    F70_ye[n] = dctrl[dd];
    dd++;
    F70_zs[n] = dctrl[dd];
    dd++;
    F70_ze[n] = dctrl[dd];
    dd++;
    }
	
	for(n=0;n<F71;++n)
    {
    F71_xs[n] = dctrl[dd];
    dd++;
    F71_xe[n] = dctrl[dd];
    dd++;
    F71_ys[n] = dctrl[dd];
    dd++;
    F71_ye[n] = dctrl[dd];
    dd++;
    F71_zs[n] = dctrl[dd];
    dd++;
    F71_ze[n] = dctrl[dd];
    dd++;
    }
	
	for(n=0;n<F72;++n)
    {
    F72_xs[n] = dctrl[dd];
    dd++;
    F72_xe[n] = dctrl[dd];
    dd++;
    F72_ys[n] = dctrl[dd];
    dd++;
    F72_ye[n] = dctrl[dd];
    dd++;
    F72_h[n] = dctrl[dd];
    dd++;
    }
	
    for(n=0;n<F369;++n)
    {
    F369_x[n] = dctrl[dd];
    dd++;
    F369_z[n] = dctrl[dd];
    dd++;
    F369_a[n] = dctrl[dd];
    dd++;
    F369_s[n] = dctrl[dd];
    dd++;
    F369_l[n] = dctrl[dd];
    dd++;
	F369_v[n] = dctrl[dd];
    dd++;
    }
    
	for(n=0;n<F370;++n)
    {
    F370_xs[n] = dctrl[dd];
    dd++;
    F370_xe[n] = dctrl[dd];
    dd++;
    F370_ys[n] = dctrl[dd];
    dd++;
    F370_ye[n] = dctrl[dd];
    dd++;
    F370_zs[n] = dctrl[dd];
    dd++;
    F370_ze[n] = dctrl[dd];
    dd++;
    }
	
	for(n=0;n<F371;++n)
    {
    F371_xs[n] = dctrl[dd];
    dd++;
    F371_xe[n] = dctrl[dd];
    dd++;
    F371_ys[n] = dctrl[dd];
    dd++;
    F371_ye[n] = dctrl[dd];
    dd++;
    F371_zs[n] = dctrl[dd];
    dd++;
    F371_ze[n] = dctrl[dd];
    dd++;
    }
    
    for(n=0;n<F374;++n)
    {
    F374_xc[n] = dctrl[dd];
    dd++;
    F374_zc[n] = dctrl[dd];
    dd++;
    F374_r[n] = dctrl[dd];
    dd++;
    }
	
    for(n=0;n<F375;++n)
    {
    F375_xc[n] = dctrl[dd];
    dd++;
    F375_zc[n] = dctrl[dd];
    dd++;
    F375_r[n] = dctrl[dd];
    dd++;
    }
    
    for(n=0;n<F378;++n)
    {
    F378_xc[n] = dctrl[dd];
    dd++;
    F378_yc[n] = dctrl[dd];
    dd++;
    F378_zc[n] = dctrl[dd];
    dd++;
    F378_r[n] = dctrl[dd];
    dd++;
    }
    
    for(n=0;n<F379;++n)
    {
    F379_xc[n] = dctrl[dd];
    dd++;
    F379_yc[n] = dctrl[dd];
    dd++;
    F379_zc[n] = dctrl[dd];
    dd++;
    F379_r[n] = dctrl[dd];
    dd++;
    }
	
	
	for(n=0;n<F390;++n)
    {
    F390_xs[n] = dctrl[dd];
    dd++;
    F390_xe[n] = dctrl[dd];
    dd++;
    F390_ys[n] = dctrl[dd];
    dd++;
    F390_ye[n] = dctrl[dd];
    dd++;
    F390_zs[n] = dctrl[dd];
    dd++;
    F390_ze[n] = dctrl[dd];
    dd++;
    }
	
	for(n=0;n<F391;++n)
    {
    F391_xs[n] = dctrl[dd];
    dd++;
    F391_xe[n] = dctrl[dd];
    dd++;
    F391_ys[n] = dctrl[dd];
    dd++;
    F391_ye[n] = dctrl[dd];
    dd++;
    F391_zs[n] = dctrl[dd];
    dd++;
    F391_ze[n] = dctrl[dd];
    dd++;
    }
    
    for(n=0;n<F394;++n)
    {
    F394_xc[n] = dctrl[dd];
    dd++;
    F394_zc[n] = dctrl[dd];
    dd++;
    F394_r[n] = dctrl[dd];
    dd++;
    }
	
    for(n=0;n<F395;++n)
    {
    F395_xc[n] = dctrl[dd];
    dd++;
    F395_zc[n] = dctrl[dd];
    dd++;
    F395_r[n] = dctrl[dd];
    dd++;
    }
    
    for(n=0;n<F398;++n)
    {
    F398_xc[n] = dctrl[dd];
    dd++;
    F398_yc[n] = dctrl[dd];
    dd++;
    F398_zc[n] = dctrl[dd];
    dd++;
    F398_r[n] = dctrl[dd];
    dd++;
    }
    
    for(n=0;n<F399;++n)
    {
    F399_xc[n] = dctrl[dd];
    dd++;
    F399_yc[n] = dctrl[dd];
    dd++;
    F399_zc[n] = dctrl[dd];
    dd++;
    F399_r[n] = dctrl[dd];
    dd++;
    }

    for(n=0;n<G60;++n)
    {
    G60_xs[n] = dctrl[dd];
    dd++;
    G60_xe[n] = dctrl[dd];
    dd++;
    G60_ys[n] = dctrl[dd];
    dd++;
    G60_ye[n] = dctrl[dd];
    dd++;
    G60_zs[n] = dctrl[dd];
    dd++;
    G60_ze[n] = dctrl[dd];
    dd++;
    }

    for(n=0;n<G61;++n)
    {
    G61_xs[n] = dctrl[dd];
    dd++;
    G61_xe[n] = dctrl[dd];
    dd++;
    G61_ys[n] = dctrl[dd];
    dd++;
    G61_ye[n] = dctrl[dd];
    dd++;
    G61_zs[n] = dctrl[dd];
    dd++;
    G61_ze[n] = dctrl[dd];
    dd++;
    }
	
	for(n=0;n<G81;++n)
    {
    G81_xs[n] = dctrl[dd];
    dd++;
    G81_xe[n] = dctrl[dd];
    dd++;
    G81_ys[n] = dctrl[dd];
    dd++;
    G81_ye[n] = dctrl[dd];
    dd++;
    G81_zs[n] = dctrl[dd];
    dd++;
    G81_ze[n] = dctrl[dd];
    dd++;
    }
	
	for(n=0;n<G95;++n)
    {
    G95_xs[n] = dctrl[dd];
    dd++;
    G95_xe[n] = dctrl[dd];
    dd++;
    G95_ys[n] = dctrl[dd];
    dd++;
    G95_ye[n] = dctrl[dd];
    dd++;
    G95_zs[n] = dctrl[dd];
    dd++;
    G95_ze[n] = dctrl[dd];
    dd++;
    }
	
	for(n=0;n<P35;++n)
    {
    P35_ts[n] = dctrl[dd];
    dd++;
    P35_te[n] = dctrl[dd];
    dd++;
	P35_dt[n] = dctrl[dd];
    dd++;
    }
	
	for(n=0;n<P50;++n)
    {
    P50_x[n] = dctrl[dd];
    dd++;
    P50_y[n] = dctrl[dd];
    dd++;
    }
	
	for(n=0;n<P51;++n)
    {
    P51_x[n] = dctrl[dd];
    dd++;
    P51_y[n] = dctrl[dd];
    dd++;
    }

    for(n=0;n<P52;++n)
    {
    P52_y[n] = dctrl[dd];
    dd++;
    }
	
	for(n=0;n<P56;++n)
    {
    P56_x[n] = dctrl[dd];
    dd++;
    }
	
	for(n=0;n<P61;++n)
    {
    P61_x[n] = dctrl[dd];
    dd++;
    P61_y[n] = dctrl[dd];
    dd++;
	P61_z[n] = dctrl[dd];
    dd++;
    }
		
	for(n=0;n<P62;++n)
    {
    P62_xs[n] = dctrl[dd];
    dd++;
    P62_xe[n] = dctrl[dd];
    dd++;
	P62_ys[n] = dctrl[dd];
    dd++;
	P62_ye[n] = dctrl[dd];
    dd++;
    P62_zs[n] = dctrl[dd];
    dd++;
	P62_ze[n] = dctrl[dd];
    dd++;
    }
	
	for(n=0;n<P67;++n)
    {
    P67_x[n] = dctrl[dd];
    dd++;
    }
	
	for(n=0;n<P81;++n)
    {
    P81_xs[n] = dctrl[dd];
    dd++;
    P81_xe[n] = dctrl[dd];
    dd++;
	P81_ys[n] = dctrl[dd];
    dd++;
	P81_ye[n] = dctrl[dd];
    dd++;
    P81_zs[n] = dctrl[dd];
    dd++;
	P81_ze[n] = dctrl[dd];
    dd++;
    }
	
	for(n=0;n<P85;++n)
    {
    P85_xs[n] = dctrl[dd];
    dd++;
    P85_xe[n] = dctrl[dd];
    dd++;
	P85_ys[n] = dctrl[dd];
    dd++;
	P85_ye[n] = dctrl[dd];
    dd++;
    P85_zs[n] = dctrl[dd];
    dd++;
	P85_ze[n] = dctrl[dd];
    dd++;
    }
	
	for(n=0;n<P121;++n)
    {
    P121_x[n] = dctrl[dd];
    dd++;
    P121_y[n] = dctrl[dd];
    dd++;
    }
	
	for(n=0;n<P123;++n)
    {
    P123_y[n] = dctrl[dd];
    dd++;
    }
	
	for(n=0;n<P124;++n)
    {
    P124_x[n] = dctrl[dd];
    dd++;
    }
	
	for(n=0;n<P125;++n)
    {
    P125_x[n] = dctrl[dd];
    dd++;
    P125_y[n] = dctrl[dd];
    dd++;
    }
    
    for(n=0;n<P230;++n)
    {
    P230_x[n] = dctrl[dd];
    dd++;
    }
    
    for(n=0;n<P240;++n)
    {
    P240_x[n] = dctrl[dd];
    dd++;
    }
	
	for(n=0;n<P351;++n)
    {
    P351_x[n] = dctrl[dd];
    dd++;
    P351_y[n] = dctrl[dd];
    dd++;
    }
	
	for(n=0;n<P352;++n)
    {
    P352_x[n] = dctrl[dd];
    dd++;
    P352_y[n] = dctrl[dd];
    dd++;
    }
    
    for(n=0;n<S73;++n)
    {
    S73_val[n]= dctrl[dd];
    dd++;
	S73_dist[n]= dctrl[dd];
    dd++;
	S73_b[n]  = dctrl[dd];
    dd++;
    S73_x[n]  = dctrl[dd];
    dd++;
	S73_y[n]  = dctrl[dd];
    dd++;
    }
	
	for(n=0;n<X110;++n)
    {
    X110_xs[n] = dctrl[dd];
    dd++;
    X110_xe[n] = dctrl[dd];
    dd++;
    X110_ys[n] = dctrl[dd];
    dd++;
    X110_ye[n] = dctrl[dd];
    dd++;
    X110_zs[n] = dctrl[dd];
    dd++;
    X110_ze[n] = dctrl[dd];
    dd++;
    }
    
    for(n=0;n<X163;++n)
    {
    X163_x1[n] = dctrl[dd];
    dd++;
    X163_y1[n] = dctrl[dd];
    dd++;
    X163_z1[n] = dctrl[dd];
    dd++;
    X163_x2[n] = dctrl[dd];
    dd++;
    X163_y2[n] = dctrl[dd];
    dd++;
    X163_z2[n] = dctrl[dd];
    dd++;
    X163_x3[n] = dctrl[dd];
    dd++;
    X163_y3[n] = dctrl[dd];
    dd++;
    X163_z3[n] = dctrl[dd];
    dd++;
    X163_x4[n] = dctrl[dd];
    dd++;
    X163_y4[n] = dctrl[dd];
    dd++;
    X163_z4[n] = dctrl[dd];
    dd++;
    X163_x5[n] = dctrl[dd];
    dd++;
    X163_y5[n] = dctrl[dd];
    dd++;
    X163_z5[n] = dctrl[dd];
    dd++;
    X163_x6[n] = dctrl[dd];
    dd++;
    X163_y6[n] = dctrl[dd];
    dd++;
    X163_z6[n] = dctrl[dd];
    dd++;
    }
    
    for(n=0;n<X164;++n)
    {
    X164_x1[n] = dctrl[dd];
    dd++;
    X164_y1[n] = dctrl[dd];
    dd++;
    X164_z1[n] = dctrl[dd];
    dd++;
    X164_x2[n] = dctrl[dd];
    dd++;
    X164_y2[n] = dctrl[dd];
    dd++;
    X164_z2[n] = dctrl[dd];
    dd++;
    X164_x3[n] = dctrl[dd];
    dd++;
    X164_y3[n] = dctrl[dd];
    dd++;
    X164_z3[n] = dctrl[dd];
    dd++;
    X164_x4[n] = dctrl[dd];
    dd++;
    X164_y4[n] = dctrl[dd];
    dd++;
    X164_z4[n] = dctrl[dd];
    dd++;
    X164_x5[n] = dctrl[dd];
    dd++;
    X164_y5[n] = dctrl[dd];
    dd++;
    X164_z5[n] = dctrl[dd];
    dd++;
    X164_x6[n] = dctrl[dd];
    dd++;
    X164_y6[n] = dctrl[dd];
    dd++;
    X164_z6[n] = dctrl[dd];
    dd++;
    X164_x7[n] = dctrl[dd];
    dd++;
    X164_y7[n] = dctrl[dd];
    dd++;
    X164_z7[n] = dctrl[dd];
    dd++;
    X164_x8[n] = dctrl[dd];
    dd++;
    X164_y8[n] = dctrl[dd];
    dd++;
    X164_z8[n] = dctrl[dd];
    dd++;
    }
    
    for(n=0;n<X311;++n)
    {
    X311_xs[n] = dctrl[dd];
    dd++;
    X311_xe[n] = dctrl[dd];
    dd++;
    X311_ys[n] = dctrl[dd];
    dd++;
    X311_ye[n] = dctrl[dd];
    dd++;
    X311_zs[n] = dctrl[dd];
    dd++;
    X311_ze[n] = dctrl[dd];
    dd++;
    X311_w[n] = dctrl[dd];
    dd++;
    X311_rho_c[n] = dctrl[dd];
    dd++;
    X311_EA[n] = dctrl[dd];
    dd++;
    X311_d[n] = dctrl[dd];
    dd++;    
    X311_l[n] = dctrl[dd];
    dd++; 
    X311_H[n] = dctrl[dd];
    dd++;
    X311_P[n] = dctrl[dd];
    dd++;
    X311_facT[n] = dctrl[dd];
    dd++;
    }
	
	//if(mpirank==1)
	//cout<<"RECV  ii: "<<ii<<" dd: "<<dd<<endl;
}
