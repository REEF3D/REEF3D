/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2022 Hans Bihs

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
Author: Hans Bihs
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
    A217 = ictrl[ii];
	ii++;
    A218 = ictrl[ii];
	ii++;
    A219 = ictrl[ii];
	ii++;
	A220 = ictrl[ii];
	ii++;
    A221 = ictrl[ii];
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
    A250 = dctrl[dd];
	dd++;
    A251 = ictrl[ii];
	ii++;
    A260 = ictrl[ii];
	ii++;
    A261 = dctrl[dd];
	dd++;
    A262 = dctrl[dd];
	dd++;
    
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
    A322 = ictrl[ii];
	ii++;
    A323 = ictrl[ii];
	ii++;
    A329 = ictrl[ii];
	ii++;
    A340 = dctrl[dd];
	dd++;
    A341 = dctrl[dd];
	dd++;
    A342 = dctrl[dd];
	dd++;
    A343 = ictrl[ii];
	ii++;
    A344 = ictrl[ii];
	ii++;
    A344_val = dctrl[dd];
	dd++;
    A345 = ictrl[ii];
	ii++;
    A345_val = dctrl[dd];
	dd++;
    A346 = dctrl[dd];
	dd++;
    A347 = ictrl[ii];
	ii++;
    A348 = ictrl[ii];
	ii++;
    A350 = ictrl[ii];
	ii++;
    A351 = ictrl[ii];
	ii++;
    A352 = ictrl[ii];
	ii++;
    A353 = ictrl[ii];
	ii++;
    A354 = dctrl[dd];
	dd++;
    A355 = dctrl[dd];
	dd++;
    A356 = dctrl[dd];
	dd++;
    A357 = ictrl[ii];
	ii++;
    A361 = ictrl[ii];
	ii++;
    A362 = ictrl[ii];
	ii++;
    A363 = ictrl[ii];
	ii++;
    A365 = dctrl[dd];
	dd++;
    A368 = ictrl[ii];
	ii++;
    A369 = dctrl[dd];
	dd++;
    
    A410 = ictrl[ii];
	ii++;
    A440 = dctrl[dd];
	dd++;
    
    A501 = ictrl[ii];
	ii++;
    A510 = ictrl[ii];
	ii++;
    A511 = ictrl[ii];
	ii++;
    A512 = ictrl[ii];
	ii++;
    A514 = ictrl[ii];
	ii++;
    A515 = ictrl[ii];
	ii++;
    A516 = ictrl[ii];
	ii++;
    A517 = ictrl[ii];
	ii++;
    A518 = ictrl[ii];
	ii++;
    A520 = ictrl[ii];
	ii++;
    A521 = ictrl[ii];
	ii++;
    A523 = dctrl[dd];
	dd++;
    A531 = dctrl[dd];
	dd++;
    A540 = ictrl[ii];
	ii++;
    A541 = dctrl[dd];
	dd++;
    A542 = dctrl[dd];
	dd++;
    A543 = ictrl[ii];
	ii++;
    A544 = dctrl[dd];
	dd++;
    A550 = ictrl[ii];
	ii++;
    A551 = ictrl[ii];
	ii++;
    A552 = ictrl[ii];
	ii++;
    A553 = ictrl[ii];
	ii++;
	
	
    B10 = ictrl[ii];
	ii++;
    B20 = ictrl[ii];
	ii++;
    B23 = ictrl[ii];
	ii++;
	B29 = dctrl[dd];
	dd++;
    B30 = ictrl[ii];
	ii++;
    B31 = dctrl[dd];
	dd++;
    B32 = ictrl[ii];
	ii++;
    B32_x = dctrl[dd];
	dd++;
    B32_y = dctrl[dd];
	dd++;
    B32_z = dctrl[dd];
	dd++;
    B33 = ictrl[ii];
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
	B70 = ictrl[ii];
	ii++;
	B71 = ictrl[ii];
	ii++;
    B75 = ictrl[ii];
	ii++;
    B76 = ictrl[ii];
	ii++;
	B77 = ictrl[ii];
	ii++;
	B81 = ictrl[ii];
	ii++;
	B81_1 = dctrl[dd];
	dd++;
	B81_2 = dctrl[dd];
	dd++;
    B81_3 = dctrl[dd];
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
    B92 = ictrl[ii];
	ii++;
    B93 = ictrl[ii];
	ii++;
    B93_1 = dctrl[dd];
	dd++;
    B93_2 = dctrl[dd];
	dd++;
    B94 = ictrl[ii];
	ii++;
    B94_wdt = dctrl[dd];
	dd++;
    B96_1 = dctrl[dd];
	dd++;
    B96_2 = dctrl[dd];
	dd++;
    B98 = ictrl[ii];
	ii++;
    B99 = ictrl[ii];
	ii++;
    B101 = ictrl[ii];
	ii++;
    B102 = dctrl[dd];
	dd++;
    B105 = ictrl[ii];
	ii++;
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
    B110 = ictrl[ii];
	ii++;
	B110_zs = dctrl[dd];
	dd++;
	B110_ze = dctrl[dd];
	dd++;
    B108 = ictrl[ii];
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
    B120 = dctrl[dd];
	dd++;
	B122 = dctrl[dd];
	dd++;
    B123 = dctrl[dd];
	dd++;
    B125 = ictrl[ii];
	ii++;
    B125_y = dctrl[dd];
	dd++;
    B127 = ictrl[ii];
	ii++;
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
    B134 = dctrl[dd];
	dd++;
    B135 = dctrl[dd];
	dd++;
    B136 = ictrl[ii];
	ii++;
    B138 = ictrl[ii];
	ii++;
    B138_1 = ictrl[ii];
	ii++;
    B138_2 = ictrl[ii];
	ii++;
    B139 = ictrl[ii];
	ii++;
    B140_1 = dctrl[dd];
	dd++;
    B140_2 = dctrl[dd];
	dd++;
    B140_3 = dctrl[dd];
	dd++;
    B160 = ictrl[ii];
	ii++;
    B170 = ictrl[ii];
	ii++;
    B180 = ictrl[ii];
	ii++;
    B181 = ictrl[ii];
	ii++;
    B181_1 = dctrl[dd];
	dd++;
    B181_2 = dctrl[dd];
	dd++;
    B181_3 = dctrl[dd];
	dd++;
    B182 = ictrl[ii];
	ii++;
    B182_1 = dctrl[dd];
	dd++;
    B182_2 = dctrl[dd];
	dd++;
    B182_3 = dctrl[dd];
	dd++;
    B183 = ictrl[ii];
	ii++;
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
    B282 = ictrl[ii];
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
    B321 = ictrl[ii];
	ii++;
    B322 = ictrl[ii];
	ii++;
    B411 = ictrl[ii];
	ii++;
    B412 = ictrl[ii];
	ii++;
    B413 = ictrl[ii];
	ii++;
    B414 = ictrl[ii];
	ii++;
    B415 = ictrl[ii];
	ii++;
    B416 = ictrl[ii];
	ii++;
    B417 = ictrl[ii];
	ii++;
    B418 = ictrl[ii];
	ii++;
    B421 = ictrl[ii];
	ii++;
    B422 = ictrl[ii];
	ii++;
    B440 = ictrl[ii];
	ii++;
    B441 = ictrl[ii];
	ii++;
    B442 = ictrl[ii];
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
    C9 = ictrl[ii];
	ii++;
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
    D30 = ictrl[ii];
	ii++;
    D31 = ictrl[ii];
	ii++;
    D32 = ictrl[ii];
	ii++;
    D33 = ictrl[ii];
	ii++;
    D37 = ictrl[ii];
	ii++;
	
    F10 = ictrl[ii];
	ii++;
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
    G2  = ictrl[ii];
	ii++;
    G3  = ictrl[ii];
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

    H1 = dctrl[dd];
	dd++;
    H2 = dctrl[dd];
	dd++;
    H3 = ictrl[ii];
	ii++;
    H4 = ictrl[ii];
	ii++;
    H4_beta1 = dctrl[dd];
	dd++;
    H4_beta2 = dctrl[dd];
	dd++;
    H9 = ictrl[ii];
	ii++;
    H10 = ictrl[ii];
	ii++;
    H15 = ictrl[ii];
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
    H61 = ictrl[ii];
	ii++;
    H61_T = dctrl[dd];
	dd++;
    H62 = ictrl[ii];
	ii++;
    H62_T = dctrl[dd];
	dd++;
    H63 = ictrl[ii];
	ii++;
    H63_T = dctrl[dd];
	dd++;
    H64 = ictrl[ii];
	ii++;
    H64_T = dctrl[dd];
	dd++;
    H65 = ictrl[ii];
	ii++;
    H65_T = dctrl[dd];
	dd++;
    H66 = ictrl[ii];
	ii++;
    H66_T = dctrl[dd];
	dd++;


    I10 = ictrl[ii];
	ii++;
    I11 = ictrl[ii];
	ii++;
    I12 = ictrl[ii];
	ii++;
    I13 = ictrl[ii];
	ii++;
    I21 = ictrl[ii];
	ii++;
	I30 = ictrl[ii];
	ii++;
	I40 = ictrl[ii];
	ii++;
	I41 = ictrl[ii];
	ii++;
    I44 = ictrl[ii];
	ii++;
    I50 = dctrl[dd];
	dd++;
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
    I233 = dctrl[dd];
	dd++;
    I240 = ictrl[ii];
	ii++;
    I241 = dctrl[dd];
	dd++;

    M10 = ictrl[ii];
	ii++;

	N10 = ictrl[ii];
    ii++;
    N11 = ictrl[ii];
    ii++;
    N40 = ictrl[ii];
    ii++;
    N41 = dctrl[dd];
    dd++;
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
    N60 = ictrl[ii];
    ii++;
    N61 = dctrl[dd];
    dd++;
    

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
	P18 = ictrl[ii];
	ii++;
	P20 = ictrl[ii];
	ii++;
    P21 = ictrl[ii];
	ii++;
    P22 = dctrl[dd];
	dd++;
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
    P43 = ictrl[ii];
	ii++;
	P43_xs = dctrl[dd];
	dd++;
    P43_xe = dctrl[dd];
	dd++;
    P43_ys = dctrl[dd];
	dd++;
    P43_ye = dctrl[dd];
	dd++;
    P44 = ictrl[ii];
	ii++;
    P45 = ictrl[ii];
	ii++;
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
    P57 = ictrl[ii];
	ii++;
    P58 = ictrl[ii];
	ii++;
	P59 = ictrl[ii];
	ii++;
	P61 = ictrl[ii];
	ii++;
	P62 = ictrl[ii];
	ii++;
    P63 = ictrl[ii];
	ii++;
    P64 = ictrl[ii];
	ii++;
	P66 = ictrl[ii];
	ii++;
	P67 = ictrl[ii];
	ii++;
    P68 = ictrl[ii];
	ii++;
    P71 = ictrl[ii];
	ii++;
    P72 = ictrl[ii];
	ii++;
    P73 = ictrl[ii];
	ii++;
    P74 = ictrl[ii];
	ii++;
    P75 = ictrl[ii];
	ii++;
	P76 = ictrl[ii];
	ii++;
    P77 = ictrl[ii];
	ii++;
    P78 = ictrl[ii];
	ii++;
    P79 = ictrl[ii];
	ii++;
    P81 = ictrl[ii];
	ii++;
    P82 = ictrl[ii];
	ii++;
	P85 = ictrl[ii];
	ii++;
	P91 = dctrl[dd];
	dd++;
    P92 = ictrl[ii];
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
	P110 = ictrl[ii];
	ii++;
    P111 = dctrl[dd];
	dd++;
    P120 = ictrl[ii];
	ii++;
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
    P184 = ictrl[ii];
	ii++;
    P185 = ictrl[ii];
	ii++;
    P190 = ictrl[ii];
	ii++;
	P191 = ictrl[ii];
	ii++;
	P192 = dctrl[dd];
	dd++;
    P194 = ictrl[ii];
	ii++;
    P195 = ictrl[ii];
	ii++;
    P230 = ictrl[ii];
	ii++;
    P240 = ictrl[ii];
	ii++;
	P351 = ictrl[ii];
	ii++;
	P352 = ictrl[ii];
	ii++;

    Q10 = ictrl[ii];
	ii++;
    Q21 = dctrl[dd];
	dd++;
    Q22 = dctrl[dd];
	dd++;
    Q23 = dctrl[dd];
	dd++;
    Q24 = ictrl[ii];
	ii++;
    Q25 = dctrl[dd];
	dd++;
    Q29 = ictrl[ii];
	ii++;
    Q31 = dctrl[dd];
	dd++;
    Q41 = dctrl[dd];
	dd++;
    Q43 = ictrl[ii];
    ii++;
    Q101 = ictrl[ii];
	ii++;
    Q110 = ictrl[ii];
	ii++;
    Q111 = ictrl[ii];
	ii++;
    Q111_x = dctrl[dd];
	dd++;
    Q112 = ictrl[ii];
	ii++;
    Q112_y = dctrl[dd];
	dd++;
    Q113 = ictrl[ii];
	ii++;
    Q113_z = dctrl[dd];
	dd++;
    Q180 = ictrl[ii];
    ii++;
    Q181 = ictrl[ii];
    ii++;
    Q182 = dctrl[dd];
	dd++;

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
	S19 = dctrl[dd];
	dd++;
    S20 = dctrl[dd];
	dd++;
	S21 = dctrl[dd];
	dd++;
    S22 = dctrl[dd];
	dd++;
    S23 = ictrl[ii];
	ii++;
    S23_val = dctrl[dd];
	dd++;
    S24 = dctrl[dd];
	dd++;
    S26_a = dctrl[dd];
	dd++;
    S26_b = dctrl[dd];
	dd++;
    S27 = ictrl[ii];
	ii++;
    S30 = dctrl[dd];
	dd++;
    S32 = ictrl[ii];
	ii++;
    S33 = ictrl[ii];
	ii++;
    S34 = ictrl[ii];
	ii++;
    S37 = ictrl[ii];
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
    S78 = ictrl[ii];
	ii++;
    S79 = ictrl[ii];
	ii++;
    S80 = ictrl[ii];
	ii++;
    S81 = dctrl[dd];
	dd++;
    S82 = dctrl[dd];
	dd++;
    S83 = ictrl[ii];
	ii++;
    S84 = ictrl[ii];
	ii++;
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
    S116 = dctrl[dd];
	dd++;
    
    T10 = ictrl[ii];
    ii++;
    T11 = ictrl[ii];
    ii++;
    T12 = ictrl[ii];
    ii++;
    T21 = ictrl[ii];
    ii++;
    T31 = dctrl[dd];
    dd++;
    T32 = dctrl[dd];
    dd++;
    T33 = ictrl[ii];
    ii++;
	T35 = dctrl[dd];
    dd++;
	T36 = ictrl[ii];
    ii++;
	T37 = dctrl[dd];
    dd++;
    T38 = dctrl[dd];
    dd++;
    T39 = ictrl[ii];
    ii++;
    T41 = ictrl[ii];
    ii++;
    T42 = dctrl[dd];
    dd++;
    T43 = dctrl[dd];
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
    W11 = ictrl[ii];
    ii++;
    W11_u = dctrl[dd];
    dd++;
    W11_v = dctrl[dd];
    dd++;
    W11_w = dctrl[dd];
    dd++;
    W12 = ictrl[ii];
    ii++;
    W12_u = dctrl[dd];
    dd++;
    W12_v = dctrl[dd];
    dd++;
    W12_w = dctrl[dd];
    dd++;
    W13 = ictrl[ii];
    ii++;
    W13_u = dctrl[dd];
    dd++;
    W13_v = dctrl[dd];
    dd++;
    W13_w = dctrl[dd];
    dd++;
    W14 = ictrl[ii];
    ii++;
    W14_u = dctrl[dd];
    dd++;
    W14_v = dctrl[dd];
    dd++;
    W14_w = dctrl[dd];
    dd++;
    W15 = ictrl[ii];
    ii++;
    W15_u = dctrl[dd];
    dd++;
    W15_v = dctrl[dd];
    dd++;
    W15_w = dctrl[dd];
    dd++;
    W16 = ictrl[ii];
    ii++;
    W16_u = dctrl[dd];
    dd++;
    W16_v = dctrl[dd];
    dd++;
    W16_w = dctrl[dd];
    dd++;
    W20 = dctrl[dd];
    dd++;
    W21 = dctrl[dd];
    dd++;
    W22 = dctrl[dd];
    dd++;
    W29_x = dctrl[dd];
    dd++;
    W29_y = dctrl[dd];
    dd++;
    W29_z = dctrl[dd];
    dd++;
	W30 = ictrl[ii];
    ii++;
	W31 = dctrl[dd];
    dd++;
    W41 = ictrl[ii];
    ii++;
    W50 = dctrl[dd];
    dd++;
    W50_air = ictrl[ii];
    ii++;
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
    X14 = ictrl[ii];
	ii++;
    X15 = ictrl[ii];
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
    X26_Cu = dctrl[dd];
	dd++;
	X26_Cv = dctrl[dd];
	dd++;
	X26_Cw = dctrl[dd];
	dd++;
    X31 = ictrl[ii];
	ii++;
	X32 = ictrl[ii];
	ii++;
	X33 = ictrl[ii];
	ii++;
    X34 = ictrl[ii];
	ii++;
    X39 = ictrl[ii];
	ii++;
    X40 = ictrl[ii];
	ii++;
	X41 = dctrl[dd];
	dd++;
	X42 = dctrl[dd];
	dd++;
	X43 = dctrl[dd];
	dd++;
	X44 = dctrl[dd];
	dd++;
    X45 = ictrl[ii];
	ii++;
    X46 = ictrl[ii];
	ii++;
    X47 = ictrl[ii];
	ii++;
    X48 = ictrl[ii];
	ii++;
    X49 = ictrl[ii];
	ii++;
    X50 = ictrl[ii];
	ii++;
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
    X181 = ictrl[ii];
	ii++;
    X181_x = dctrl[dd];
	dd++;
    X181_y = dctrl[dd];
	dd++;
    X181_z = dctrl[dd];
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
    X184 = dctrl[dd];
	dd++;
    X205 = ictrl[ii];
	ii++;
    X206 = ictrl[ii];
	ii++;
	X206_ts = dctrl[dd];
	dd++;
    X206_te = dctrl[dd];
	dd++;
    X207 = ictrl[ii];
	ii++;
	X207_ts = dctrl[dd];
	dd++;
    X207_te = dctrl[dd];
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
    X312 = ictrl[ii];
	ii++;    
    X313 = ictrl[ii];
	ii++;    
    X314 = ictrl[ii];
	ii++;    
    X315 = ictrl[ii];
	ii++;    
    X320 = ictrl[ii];
	ii++;
    X321 = ictrl[ii];
	ii++;
    X323_m = dctrl[dd];
	dd++;
    X323_d = dctrl[dd];
	dd++;
    X323_l = dctrl[dd];
	dd++;
	X324 = ictrl[ii];
	ii++;
    X325_dt = dctrl[dd];
	dd++;
    X325_relX = dctrl[dd];
	dd++;
    X325_relY = dctrl[dd];
	dd++;
    X325_relZ = dctrl[dd];
	dd++;
	X400 = ictrl[ii];
	ii++;
	X401_p0 = dctrl[dd];
    dd++;
    X401_cl = dctrl[dd]; 
    dd++;
    X401_cb = dctrl[dd]; 
    dd++;
    X401_a = dctrl[dd]; 
    dd++;
	

    Y1 = ictrl[ii];
	ii++;
    Y2 = ictrl[ii];
	ii++;
    Y3 = ictrl[ii];
	ii++;
    Y4 = ictrl[ii];
	ii++;
    Y5 = ictrl[ii];
	ii++;
    Y40 = ictrl[ii];
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
    
	Z10 = ictrl[ii];
	ii++;
    Z11 = ictrl[ii];
    ii++;
    Z12_cdx = dctrl[dd];
	dd++;
    Z12_cdy = dctrl[dd];
	dd++;
    Z12_cdz = dctrl[dd];
	dd++;
    Z12_ckx = dctrl[dd];
	dd++;
    Z12_cky = dctrl[dd];
	dd++;
    Z12_ckz = dctrl[dd];
	dd++;

// --------------------------	
	
	
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

    if(B282>0)
	{	
	Darray(B282_xs,B282);
	Darray(B282_xe,B282);
	Darray(B282_ys,B282);
	Darray(B282_ye,B282);
	Darray(B282_zs,B282);
	Darray(B282_ze,B282);
    Darray(B282_n,B282);
    Darray(B282_d50,B282);
	Darray(B282_alpha,B282);
	Darray(B282_beta,B282);
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

    if(B321>0)
	{	
	Darray(B321_xs,B321);
	Darray(B321_xe,B321);
	Darray(B321_ys,B321);
	Darray(B321_ye,B321);
	Darray(B321_zs,B321);
	Darray(B321_ze,B321);
    Darray(B321_N,B321);
    Darray(B321_D,B321);
	Darray(B321_Cd,B321);
	}

    if(B322>0)
	{	
	Darray(B322_xs,B322);
	Darray(B322_xe,B322);
	Darray(B322_ys,B322);
	Darray(B322_ye,B322);
	Darray(B322_zs,B322);
	Darray(B322_ze,B322);
    Darray(B322_N,B322);
    Darray(B322_D,B322);
	Darray(B322_Cd,B322);
	}
    
    if(B411>0)
    {
    Iarray(B411_ID,B411);
    Darray(B411_Q,B411);
    }
    
    if(B412>0)
    {
    Iarray(B412_ID,B412);
    Darray(B412_pressBC,B412);
    }
    
    if(B413>0)
    {
    Iarray(B413_ID,B413);
    Darray(B413_h,B413);
    }
    
    if(B414>0)
    {
    Iarray(B414_ID,B414);
    Darray(B414_Uio,B414);
    }
    
    if(B415>0)
    {
    Iarray(B415_ID,B415);
    Darray(B415_U,B415);
    Darray(B415_V,B415);
    Darray(B415_W,B415);
    }
    
    if(B416>0)
    {
    Iarray(B416_ID,B416);
    Darray(B416_alpha,B416);
    }
    
    if(B417>0)
    {
    Iarray(B417_ID,B417);
    Darray(B417_Nx,B417);
    Darray(B417_Ny,B417);
    Darray(B417_Nz,B417);
    }
    
    if(B418>0)
    {
    Iarray(B418_ID,B418);
    Iarray(B418_pio,B418);
    }
    
    if(B421>0)
    {
    Iarray(B421_ID,B421);
    Iarray(B421_Q,B421);
    }
    
    if(B422>0)
    {
    Iarray(B422_ID,B422);
    Iarray(B422_FSF,B422);
    }
    
    if(B440>0)
    {
    Iarray(B440_ID,B440);
    Iarray(B440_face,B440);
    Darray(B440_xs,B440);
    Darray(B440_xe,B440);
    Darray(B440_ys,B440);
    Darray(B440_ye,B440);
    }
    
    if(B441>0)
    {
    Iarray(B441_ID,B441);
    Iarray(B441_face,B441);
    Darray(B441_xs,B441);
    Darray(B441_xe,B441);
    Darray(B441_ys,B441);
    Darray(B441_ye,B441);
    Darray(B441_zs,B441);
    Darray(B441_ze,B441);
    }
    
    if(B442>0)
    {
    Iarray(B442_ID,B442);
    Iarray(B442_face,B442);
    Darray(B442_xm,B442);
    Darray(B442_ym,B442);
    Darray(B442_zm,B442);
    Darray(B442_r,B442);
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

    if(P58>0)
	{
    Darray(P58_x,P58);  
	Darray(P58_y,P58); 
	Darray(P58_T,P58);  
	}
	
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
    
    if(P63>0)
	{
    Darray(P63_x,P63);  
	Darray(P63_y,P63); 
	}

    if(P64>0)
	{
    Darray(P64_x,P64);  
	Darray(P64_y,P64); 
	Darray(P64_z,P64);  
	}
	
	if(P67>0)
	Darray(P67_x,P67); 

    if(P68>0)
    {
	Darray(P68_x,P68);
    Darray(P68_zs,P68);
    Darray(P68_ze,P68);
    }
	
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
	Darray(P85_x,P85); 
	Darray(P85_y,P85); 
	Darray(P85_r,P85);
	Darray(P85_cd,P85);
	Darray(P85_cm,P85);
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
    
    if(P184>0)
	{
    Iarray(P184_its,P184);  
	Iarray(P184_ite,P184);  
	Iarray(P184_dit,P184);  
	}
    
    if(P185>0)
	{
    Darray(P185_ts,P185);  
	Darray(P185_te,P185);  
	Darray(P185_dt,P185);  
	}

    if(P194>0)
	{
    Iarray(P194_its,P194);  
	Iarray(P194_ite,P194);  
	Iarray(P194_dit,P194);  
	}
    
    if(P195>0)
	{
    Darray(P195_ts,P195);  
	Darray(P195_te,P195);  
	Darray(P195_dt,P195);  
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
    
    if(Q110>0)
	{
	Darray(Q110_xs,Q110);  
	Darray(Q110_xe,Q110);  
	
	Darray(Q110_ys,Q110);  
	Darray(Q110_ye,Q110);  
	
	Darray(Q110_zs,Q110);  
	Darray(Q110_ze,Q110);  
	}

    if(S73>0)
	{
	Darray(S73_val,S73);
	Darray(S73_dist,S73);
	Darray(S73_b,S73);
	Darray(S73_x,S73);
	Darray(S73_y,S73);
	}
    
    if(W41>0)
    {
    Darray(W41_xc,W41);
    Darray(W41_yc,W41);
    Darray(W41_zs,W41);
    Darray(W41_ze,W41);
    Darray(W41_vel,W41);
    Darray(W41_beta,W41);
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
		
        Darray(X314_T,X311);  
        Darray(X315_t,X311);  
	}

    if(X312>0)
	{
		Darray(X311_xs,X312);
		Darray(X311_xe,X312);
		Darray(X311_ys,X312);
		Darray(X311_ye,X312);
		Darray(X311_zs,X312);
		Darray(X311_ze,X312);
		Darray(X312_k,X312);
		Darray(X312_T0,X312);
		
        Darray(X314_T,X312);  
        Darray(X315_t,X312);  
	}

    if(X320>0)
	{
		Iarray(X320_type,X320);
    }

    if(X321>0)
	{
		Darray(X321_Sn,X321);
		Darray(X321_d,X321);
        Darray(X321_lambda,X321);
		Darray(X321_dk,X321);
		Darray(X321_rho,X321);
		Darray(X321_nd,X321);
		Darray(X321_nl,X321);

        Darray(X322_D,X321);
		Darray(X322_L,X321);
        Darray(X322_x0,X321);
		Darray(X322_y0,X321);
		Darray(X322_z0,X321);
		Darray(X322_phi,X321);
		Darray(X322_theta,X321);
		Darray(X322_psi,X321);
	}

    if(X324>0)
	{
    Darray(X324_x,X324);
	Darray(X324_y,X324);
	Darray(X324_z,X324);
	}
    
    if(Z11>0)
	{
		Darray(Z11_x,Z11);  
		Darray(Z11_y,Z11);  
		Darray(Z11_z,Z11);  
		Darray(Z11_l,Z11);  
		Darray(Z11_w,Z11);  
		Darray(Z11_t,Z11);  
		Darray(Z11_rho,Z11);  
		Darray(Z11_e,Z11);  
		Darray(Z11_ix,Z11);  
		Darray(Z11_iy,Z11);  
		Darray(Z11_iz,Z11);  
		Darray(Z11_nu,Z11);  
		Darray(Z11_n,Z11);  
	}

// --------------------------

	
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

    for(n=0;n<B282;++n)
    {
	B282_xs[n] = dctrl[dd];
    dd++;
    B282_xe[n] = dctrl[dd];
    dd++;
    B282_ys[n] = dctrl[dd];
    dd++;
    B282_ye[n] = dctrl[dd];
    dd++;
    B282_zs[n] = dctrl[dd];
    dd++;
    B282_ze[n] = dctrl[dd];
    dd++;
    B282_n[n]  = dctrl[dd];
    dd++;
    B282_d50[n]= dctrl[dd];
    dd++;
	B282_alpha[n]= dctrl[dd];
    dd++;
	B282_beta[n]= dctrl[dd];
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

    for(n=0;n<B321;++n)
    {
	B321_xs[n] = dctrl[dd];
    dd++;
    B321_xe[n] = dctrl[dd];
    dd++;
    B321_ys[n] = dctrl[dd];
    dd++;
    B321_ye[n] = dctrl[dd];
    dd++;
    B321_zs[n] = dctrl[dd];
    dd++;
    B321_ze[n] = dctrl[dd];
    dd++;
    B321_N[n]  = dctrl[dd];
    dd++;
    B321_D[n]  = dctrl[dd];
    dd++;
	B321_Cd[n] = dctrl[dd];
    dd++;
    }

    for(n=0;n<B322;++n)
    {
	B322_xs[n] = dctrl[dd];
    dd++;
    B322_xe[n] = dctrl[dd];
    dd++;
    B322_ys[n] = dctrl[dd];
    dd++;
    B322_ye[n] = dctrl[dd];
    dd++;
    B322_zs[n] = dctrl[dd];
    dd++;
    B322_ze[n] = dctrl[dd];
    dd++;
    B322_N[n]  = dctrl[dd];
    dd++;
    B322_D[n]  = dctrl[dd];
    dd++;
	B322_Cd[n] = dctrl[dd];
    dd++;
    }
    
    for(n=0;n<B411;++n)
    {
    B411_ID[n] = ictrl[ii];
    ii++;
    B411_Q[n] = dctrl[dd];
    dd++;
    }
    
    for(n=0;n<B412;++n)
    {
    B412_ID[n] = ictrl[ii];
    ii++;
    B412_pressBC[n] = dctrl[dd];
    dd++;
    }
    
    for(n=0;n<B413;++n)
    {
    B413_ID[n] = ictrl[ii];
    ii++;
    B413_h[n] = dctrl[dd];
    dd++;
    }
    
    for(n=0;n<B414;++n)
    {
    B414_ID[n] = ictrl[ii];
    ii++;
    B414_Uio[n] = dctrl[dd];
    dd++;
    }
    
    for(n=0;n<B415;++n)
    {
    B415_ID[n] = ictrl[ii];
    ii++;
    B415_U[n] = dctrl[dd];
    dd++;
    B415_V[n] = dctrl[dd];
    dd++;
    B415_W[n] = dctrl[dd];
    dd++;
    }
    
    for(n=0;n<B416;++n)
    {
    B416_ID[n] = ictrl[ii];
    ii++;
    B416_alpha[n] = dctrl[dd];
    dd++;
    }
    
    for(n=0;n<B417;++n)
    {
    B417_ID[n] = ictrl[ii];
    ii++;
    B417_Nx[n] = dctrl[dd];
    dd++;
    B417_Ny[n] = dctrl[dd];
    dd++;
    B417_Nz[n] = dctrl[dd];
    dd++;
    }
    
    for(n=0;n<B418;++n)
    {
    B418_ID[n] = ictrl[ii];
    ii++;
    B418_pio[n] = ictrl[ii];
    ii++;
    }
    
    for(n=0;n<B421;++n)
    {
    B421_ID[n] = ictrl[ii];
    ii++;
    B421_Q[n] = ictrl[ii];
    ii++;
    }
	
    for(n=0;n<B422;++n)
    {
    B422_ID[n] = ictrl[ii];
    ii++;
    B422_FSF[n] = ictrl[ii];
    ii++;
    }
	
	for(n=0;n<B440;++n)
    {
    B440_ID[n]  = ictrl[ii];
    ii++;
    B440_face[n]  = ictrl[ii];
    ii++;
	B440_xs[n] = dctrl[dd];
    dd++;
    B440_xe[n] = dctrl[dd];
    dd++;
    B440_ys[n] = dctrl[dd];
    dd++;
    B440_ye[n] = dctrl[dd];
    dd++;
    }
    
    for(n=0;n<B441;++n)
    {
    B441_ID[n]  = ictrl[ii];
    ii++;
    B441_face[n]  = ictrl[ii];
    ii++;
	B441_xs[n] = dctrl[dd];
    dd++;
    B441_xe[n] = dctrl[dd];
    dd++;
    B441_ys[n] = dctrl[dd];
    dd++;
    B441_ye[n] = dctrl[dd];
    dd++;
    B441_zs[n] = dctrl[dd];
    dd++;
    B441_ze[n] = dctrl[dd];
    dd++;
    }
    
    for(n=0;n<B442;++n)
    {
    B442_ID[n]  = ictrl[ii];
    ii++;
    B442_face[n]  = ictrl[ii];
    ii++;
	B442_xm[n] = dctrl[dd];
    dd++;
    B442_ym[n] = dctrl[dd];
    dd++;
    B442_zm[n] = dctrl[dd];
    dd++;
    B442_r[n] = dctrl[dd];
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

    for(n=0;n<P58;++n)
    {
    P58_x[n] = dctrl[dd];
    dd++;
    P58_y[n] = dctrl[dd];
    dd++;
	P58_T[n] = dctrl[dd];
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
    
    for(n=0;n<P63;++n)
    {
    P63_x[n] = dctrl[dd];
    dd++;
    P63_y[n] = dctrl[dd];
    dd++;
    }

    for(n=0;n<P64;++n)
    {
    P64_x[n] = dctrl[dd];
    dd++;
    P64_y[n] = dctrl[dd];
    dd++;
	P64_z[n] = dctrl[dd];
    dd++;
    }
	
	for(n=0;n<P67;++n)
    {
    P67_x[n] = dctrl[dd];
    dd++;
    }

    for(n=0;n<P68;++n)
    {
    P68_x[n] = dctrl[dd];
    dd++;
    P68_zs[n] = dctrl[dd];
    dd++;
    P68_ze[n] = dctrl[dd];
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
    P85_x[n] = dctrl[dd];
    dd++;
    P85_y[n] = dctrl[dd];
    dd++;
	 P85_r[n] = dctrl[dd];
    dd++;
	 P85_cd[n] = dctrl[dd];
    dd++;
	 P85_cm[n] = dctrl[dd];
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
    
    for(n=0;n<P184;++n)
    {
    P184_its[n] = ictrl[ii];
    ii++;
    P184_ite[n] = ictrl[ii];
    ii++;
	P184_dit[n] = ictrl[ii];
    ii++;
    }
    
    for(n=0;n<P185;++n)
    {
    P185_ts[n] = dctrl[dd];
    dd++;
    P185_te[n] = dctrl[dd];
    dd++;
	P185_dt[n] = dctrl[dd];
    dd++;
    }

    for(n=0;n<P194;++n)
    {
    P194_its[n] = ictrl[ii];
    ii++;
    P194_ite[n] = ictrl[ii];
    ii++;
	P194_dit[n] = ictrl[ii];
    ii++;
    }
    
    for(n=0;n<P195;++n)
    {
    P195_ts[n] = dctrl[dd];
    dd++;
    P195_te[n] = dctrl[dd];
    dd++;
	P195_dt[n] = dctrl[dd];
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

    for(n=0;n<Q110;++n)
    {
    Q110_xs[n] = dctrl[dd];
    dd++;
    Q110_xe[n] = dctrl[dd];
    dd++;
    Q110_ys[n] = dctrl[dd];
    dd++;
    Q110_ye[n] = dctrl[dd];
    dd++;
    Q110_zs[n] = dctrl[dd];
    dd++;
    Q110_ze[n] = dctrl[dd];
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
    
    for(n=0;n<W41;++n)
    {
    W41_xc[n] = dctrl[dd];
    dd++;
    W41_yc[n] = dctrl[dd];
    dd++;
    W41_zs[n] = dctrl[dd];
    dd++;
    W41_ze[n] = dctrl[dd];
    dd++;
    W41_vel[n] = dctrl[dd];
    dd++;
    W41_beta[n] = dctrl[dd];
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

    for(n=0;n<X312;++n)
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
    X312_k[n] = dctrl[dd];
    dd++;
    X312_T0[n] = dctrl[dd];
    dd++;
    }

    if (X314 > 0)
    {
        for(n=0;n<X311;++n)
        {
            X314_T[n] = dctrl[dd]; 
            dd++;
        }
        for(n=0;n<X312;++n)
        {
            X314_T[n] = dctrl[dd]; 
            dd++;
        }
    }
    if (X315 > 0)
    {
        for(n=0;n<X311;++n)
        {
            X315_t[n] = dctrl[dd];
            dd++;
        }
        for(n=0;n<X312;++n)
        {
            X315_t[n] = dctrl[dd]; 
            dd++;
        }
    }

    for(n=0;n<X321;++n)
    {
    X321_Sn[n] = dctrl[dd];
    dd++;
    X321_d[n] = dctrl[dd];
    dd++;
    X321_lambda[n] = dctrl[dd];
    dd++;
    X321_dk[n] = dctrl[dd];
    dd++;
    X321_rho[n] = dctrl[dd];
    dd++;
    X321_nd[n] = dctrl[dd];
    dd++;
    X321_nl[n] = dctrl[dd];
    dd++;

    X322_D[n] = dctrl[dd];
    dd++;
    X322_L[n] = dctrl[dd];
    dd++;
    X322_x0[n] = dctrl[dd];
    dd++;
    X322_y0[n] = dctrl[dd];
    dd++;
    X322_z0[n] = dctrl[dd];
    dd++;
    X322_phi[n] = dctrl[dd];
    dd++;
    X322_theta[n] = dctrl[dd];
    dd++;
    X322_psi[n] = dctrl[dd];
    dd++;
    }

    for(n=0;n<X324;++n)
    {
    X324_x[n] = dctrl[dd];
    dd++;
    X324_y[n] = dctrl[dd];
    dd++;
	X324_z[n] = dctrl[dd];
    dd++;
    }

    for(n=0;n<Z11;++n)
    {
    Z11_x[n] = dctrl[dd];
    dd++;
    Z11_y[n] = dctrl[dd];
    dd++;
    Z11_z[n] = dctrl[dd];
    dd++;
    Z11_l[n] = dctrl[dd];
    dd++;
    Z11_w[n] = dctrl[dd];
    dd++;
    Z11_t[n] = dctrl[dd];
    dd++;
    Z11_rho[n] = dctrl[dd];
    dd++;
    Z11_e[n] = dctrl[dd];
    dd++;
    Z11_ix[n] = dctrl[dd];
    dd++;
    Z11_iy[n] = dctrl[dd];
    dd++;
    Z11_iz[n] = dctrl[dd];
    dd++;
    Z11_nu[n] = dctrl[dd];
    dd++;
    Z11_n[n] = dctrl[dd];
    dd++;
    }
}
