// Unit1.cpp
//---------------------------------------------------------------------------
#include "Unit2.h"
#include <iostream.h>
#include <math.h>
#include <Windows.h>
#include <iomanip.h>
#include <time.h>
#include <cstring.h>
#include <cstddef.h>
#include <fstream>
#include <vcl.h>
#pragma hdrstop

#include "Unit1.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TForm1 *Form1;
double t, ny, ny_1;

double Pi = 3.14159;

double Ro0 = 1.225;

double Rz = 6371000;

double Pi_0 = 3.986005E+14;

double Hm = 7110;

double  h = 1;

double g = 9.80665;

double t_k;

double m_0;

double m_k;

double P_ud_p;

double Cx0;

double Cya;

double S;

double S_a;

double V, r, L, Vx, Vy;
double rx, ry = Rz;
double alpha, teta, teta_v, eta;
double  w1, peregruz_x, peregruz_y, P_not_cos_sin, q, H;
double P, Rx, Ry, wd1;
double a1, a2, a3, a4;
double Kx_P1, Kx_Cx1, Kx_Cy1, Kx_G1;
double Ky_P1, Ky_Cx1, Ky_Cy1, Ky_G1;
double Kx_P2, Kx_Cx2, Kx_Cy2, Kx_G2;
double Ky_P2, Ky_Cx2, Ky_Cy2, Ky_G2;
double Kx_P3, Kx_Cx3, Kx_Cy3, Kx_G3;
double Ky_P3, Ky_Cx3, Ky_Cy3, Ky_G3;
double Kx_P4, Kx_Cx4, Kx_Cy4, Kx_G4;
double Ky_P4, Ky_Cx4, Ky_Cy4, Ky_G4;
double  k1_rx, k1_ry, k1_Vx, k1_Vy;
double  k2_rx, k2_ry, k2_Vx, k2_Vy;
double  k3_rx, k3_ry, k3_Vx, k3_Vy;
double  k4_rx, k4_ry, k4_Vx, k4_Vy;
double r_old[2], r_new[2];
double del_r, del_rx, del_ry, S_old, del_r_old;
TLineSeries *NewSeries1;
TLineSeries *NewSeries2;
TLineSeries *NewSeries3;
TLineSeries *NewSeries4;
TLineSeries *NewSeries5;
TLineSeries *NewSeries6;
TLineSeries *NewSeries7;
TLineSeries *NewSeries8;
TLineSeries *NewSeries9;
TLineSeries *NewSeries10;
TLineSeries *NewSeries11;
//---------------------------------------------------------------------------
__fastcall TForm1::TForm1(TComponent* Owner)
        : TForm(Owner)
{
}
//---------------------------------------------------------------------------

void __fastcall TForm1::Button1Click(TObject *Sender)
{
   ofstream output("–езультаты моделировани€.txt");
   output << "t" << setw(12) << "w1" << setw(16) << "V" << setw(15) << "q" << setw(15)
        << "peregruz_x" << setw(14) << "peregruz_y" << setw(14) << "ny" << setw(15) << setw(14)
        << "L" << setw(15) << "H" << setw(15) << "eta" << setw(14) << setw(16) << "teta" << setw(16)
        << "alpha" << endl;
   srand(time(NULL));
   V = 0, r = 0, L = 0, Vx = 0, Vy = 0;
   rx = 0, ry = Rz;
   alpha = 0, teta = 0, teta_v = 0, eta = 0;
   w1 = 0, peregruz_x = 0, peregruz_y = 0, P_not_cos_sin = 0, q = 0, H = 0;
   P = 0, Rx = 0, Ry = 0, wd1 = 0;
   a1 = 0, a2 = 0, a3 = 0, a4 = 0;
   del_r = 0, del_rx = 0, del_ry = 0, del_r_old = 0;
   Kx_P1 = 0, Kx_Cx1 = 0, Kx_Cy1 = 0, Kx_G1 = 0;
   Ky_P1 = 0, Ky_Cx1 = 0, Ky_Cy1 = 0, Ky_G1 = 0;
   Kx_P2 = 0, Kx_Cx2 = 0, Kx_Cy2 = 0, Kx_G2 = 0;
   Ky_P2 = 0, Ky_Cx2 = 0, Ky_Cy2 = 0, Ky_G2 = 0;
   Kx_P3 = 0, Kx_Cx3 = 0, Kx_Cy3 = 0, Kx_G3 = 0;
   Ky_P3 = 0, Ky_Cx3 = 0, Ky_Cy3 = 0, Ky_G3 = 0;
   Kx_P4 = 0, Kx_Cx4 = 0, Kx_Cy4 = 0, Kx_G4 = 0;
   Ky_P4 = 0, Ky_Cx4 = 0, Ky_Cy4 = 0, Ky_G4 = 0;
   k1_rx = 0, k1_ry = 0, k1_Vx = 0, k1_Vy = 0;
   k2_rx = 0, k2_ry = 0, k2_Vx = 0, k2_Vy = 0;
   k3_rx = 0, k3_ry = 0, k3_Vx = 0, k3_Vy = 0;
   k4_rx = 0, k4_ry = 0, k4_Vx = 0, k4_Vy = 0;
   TColor OColor[18] = {RGB(28, 243, 17), RGB(0, 0, 0), RGB(8, 20, 252), RGB(181, 74, 151), RGB(52, 160, 208),
                        RGB(49, 211, 111), RGB(232, 149, 28), RGB(255, 5, 184), RGB(255, 5, 5), RGB(171, 167, 89),
                        RGB(148, 112, 112), RGB(99, 111, 161), RGB(92, 168, 106), RGB(159, 101, 156), RGB(130, 130, 130),
                        RGB(87, 120, 173), RGB(173, 87, 87), RGB(150, 110, 141)};
   int ParamColor[11] = {random(18), random(18), random(18), random(18),
                        random(18), random(18), random(18), random(18),
                        random(18), random(18), random(18)};
   if(CheckBox1 -> Checked)
   {
       NewSeries1 = new TLineSeries(Chart1);
       Chart1 -> AddSeries(NewSeries1);
       Label11 -> Visible = true;
       Label12 -> Visible = true;
       Label11 -> Caption = "\tH, м";
       Label12 -> Caption = " t, с";
   }
   if(CheckBox2 -> Checked)
   {
       NewSeries2 = new TLineSeries(Chart1);
       Chart1 -> AddSeries(NewSeries2);
       Label11 -> Visible = true;
       Label12 -> Visible = true;
       Label11 -> Caption = "\tL, м";
       Label12 -> Caption = " t, с";
   }
   if(CheckBox3 -> Checked)
   {
       NewSeries3 = new TLineSeries(Chart1);
       Chart1 -> AddSeries(NewSeries3);
       Label11 -> Visible = true;
       Label12 -> Visible = true;
       Label11 -> Caption = " alpha, рад";
       Label12 -> Caption = " t, с";
   }
   if(CheckBox4 -> Checked)
   {
       NewSeries4 = new TLineSeries(Chart1);
       Chart1 -> AddSeries(NewSeries4);
       Label11 -> Visible = true;
       Label12 -> Visible = true;
       Label11 -> Caption = "    q, кг/(м * c^2)";
       Label12 -> Caption = " t, с";
   }
   if(CheckBox5 -> Checked)
   {
       NewSeries5 = new TLineSeries(Chart1);
       Chart1 -> AddSeries(NewSeries5);
       Label11 -> Visible = true;
       Label12 -> Visible = true;
       Label11 -> Caption = "\tV, м/с";
       Label12 -> Caption = " t, с";
   }
   if(CheckBox6 -> Checked)
   {
       NewSeries6 = new TLineSeries(Chart1);
       Chart1 -> AddSeries(NewSeries6);
       Label11 -> Visible = true;
       Label12 -> Visible = true;
       Label11 -> Caption = "peregruz_x, м/с^2";
       Label12 -> Caption = " t, с";
   }
   if(CheckBox7 -> Checked)
   {
       NewSeries7 = new TLineSeries(Chart1);
       Chart1 -> AddSeries(NewSeries7);
       Label11 -> Visible = true;
       Label12 -> Visible = true;
       Label11 -> Caption = "peregruz_y, м/с^2 ";
       Label12 -> Caption = " t, с";
   }
   if(CheckBox8 -> Checked)
   {
       NewSeries8 = new TLineSeries(Chart1);
       Chart1 -> AddSeries(NewSeries8);
       Label11 -> Visible = true;
       Label12 -> Visible = true;
       Label11 -> Caption = "\tw1, м/с";
       Label12 -> Caption = " t, с";
   }
   if(CheckBox9 -> Checked)
   {
       NewSeries9 = new TLineSeries(Chart1);
       Chart1 -> AddSeries(NewSeries9);
       Label11 -> Visible = true;
       Label12 -> Visible = true;
       Label11 -> Caption = "   ny, рад";
       Label12 -> Caption = " t, с";
   }
   if(CheckBox10 -> Checked)
   {
       NewSeries10 = new TLineSeries(Chart1);
       Chart1 -> AddSeries(NewSeries10);
       Label11 -> Visible = true;
       Label12 -> Visible = true;
       Label11 -> Caption = "   eta, рад";
       Label12 -> Caption = " t, с";
   }
   if(CheckBox11 -> Checked)
   {
       NewSeries11 = new TLineSeries(Chart1);
       Chart1 -> AddSeries(NewSeries11);
       Label11 -> Visible = true;
       Label12 -> Visible = true;
       Label11 -> Caption = "    teta, рад";
       Label12 -> Caption = " t, с";
   }
   ny = StrToFloat (Edit1 -> Text)* Pi/180;
   ny_1 = StrToFloat (Edit2 -> Text)* Pi/180;
   m_0 = StrToFloat (Edit3 -> Text);
   m_k = StrToFloat (Edit4 -> Text);
   t_k = StrToFloat (Edit5 -> Text);
   P_ud_p = StrToFloat (Edit6 -> Text);
   Cx0 = StrToFloat (Edit7 -> Text);
   Cya = StrToFloat (Edit8 -> Text);
   S = StrToFloat (Edit9 -> Text);
   S_a = StrToFloat (Edit10 -> Text);
    for (int t = 0; t <= 120; t++)
    {
    if (t < 8)
       ny = ny;
    if (t >= 8 && t <= 68)
       ny -= (ny - ny_1) / 60;
    if (t > 68)
       ny = ny;
    if (Vx == 0)
       teta_v = 1.5704;
    else
       teta_v = atan(Vy / Vx);
        if (H < 100000)
        {
        P = (P_ud_p * m_k / t_k) * g;
        Rx = Cx0 * S * q;
        Ry = Cya * alpha * S * q;
        wd1 = (P + sqrt(pow(Rx, 2) + pow(Ry, 2))) / m(t, m_0, m_k, t_k);

        a1 = sqrt(pow(rx, 2) + pow(ry, 2));
        Kx_P1 = ((P_ud_p * m_k / t_k) * g) * cos(ny);
        Kx_Cx1 = Ro(H) * Cx0 * S * Vx * sqrt(pow(Vx, 2) + pow(Vy, 2)) / 2;
        Kx_Cy1 = Ro(H) * Cya * alpha * S * Vy * sqrt(pow(Vx, 2) + pow(Vy, 2)) / 2;
        Kx_G1 = (Pi_0 * rx) / pow(a1, 3);
        Ky_P1 = ((P_ud_p * m_k / t_k) * g) * sin(ny);
        Ky_Cx1 = Ro(H) * Cx0 * S * Vy * sqrt(pow(Vx, 2) + pow(Vy, 2)) / 2;
        Ky_Cy1 = Ro(H) * Cya * alpha * S * Vx * sqrt(pow(Vx, 2) + pow(Vy, 2)) / 2;
        Ky_G1 = (Pi_0 * ry) / pow(a1, 3);

        k1_rx = Vx;
        k1_ry = Vy;
        k1_Vx = (1 / m(t, m_0, m_k, t_k)) * (Kx_P1 - Kx_Cx1 - Kx_Cy1) - Kx_G1;
        k1_Vy = (1 / m(t, m_0, m_k, t_k)) * (Ky_P1 - Ky_Cx1 + Ky_Cy1) - Ky_G1;

        a2 = sqrt(pow(rx + h * k1_rx / 2, 2) + pow(ry + h * k1_ry / 2, 2));
        Kx_P2 = ((P_ud_p * m_k / t_k) * g) * cos(ny);
        Kx_Cx2 = Ro(H) * Cx0 * S * (Vx + h * k1_Vx / 2) * sqrt(pow((Vx + h * k1_Vx / 2), 2) + pow((Vy + h * k1_Vy / 2), 2)) / 2;
        Kx_Cy2 = Ro(H) * Cya * alpha * S * (Vy + h * k1_Vy / 2) * sqrt(pow((Vx + h * k1_Vx / 2), 2) + pow((Vy + h * k1_Vy / 2), 2)) / 2;
        Kx_G2 = (Pi_0 * (rx + h * k1_rx / 2)) / pow(a2, 3);
        Ky_P2 = ((P_ud_p * m_k / t_k) * g) * sin(ny);
        Ky_Cx2 = Ro(H) * Cx0 * S * (Vy + h * k1_Vy / 2) * sqrt(pow((Vx + h * k1_Vx / 2), 2) + pow((Vy + h * k1_Vy / 2), 2)) / 2;
        Ky_Cy2 = Ro(H) * Cya * alpha * S * (Vx + h * k1_Vx / 2) * sqrt(pow((Vx + h * k1_Vx / 2), 2) + pow((Vy + h * k1_Vy / 2), 2)) / 2;
        Ky_G2 = (Pi_0 * (ry + h * k1_ry)) / pow(a2, 3);

        k2_rx = Vx + h * k1_Vx / 2;
        k2_ry = Vy + h * k1_Vy / 2;
        k2_Vx = (1 / m(t, m_0, m_k, t_k)) * (Kx_P2 - Kx_Cx2 - Kx_Cy2) - Kx_G2;
        k2_Vy = (1 / m(t, m_0, m_k, t_k)) * (Ky_P2 - Ky_Cx2 + Ky_Cy2) - Ky_G2;

        a3 = sqrt(pow(rx + h * k2_rx / 2, 2) + pow(ry + h * k2_ry / 2, 2));
        Kx_P3 = ((P_ud_p * m_k / t_k) * g) * cos(ny);
        Kx_Cx3 = Ro(H) * Cx0 * S * (Vx + h * k2_Vx / 2) * sqrt(pow((Vx + h * k2_Vx / 2), 2) + pow((Vy + h * k2_Vy / 2), 2)) / 2;
        Kx_Cy3 = Ro(H) * Cya * alpha * S * (Vy + h * k2_Vy / 2) * sqrt(pow((Vx + h * k2_Vx / 2), 2) + pow((Vy + h * k2_Vy / 2), 2)) / 2;
        Kx_G3 = (Pi_0 * (rx + h * k2_rx / 2)) / pow(a3, 3);
        Ky_P3 = ((P_ud_p * m_k / t_k) * g) * sin(ny);
        Ky_Cx3 = Ro(H) * Cx0 * S * (Vy + h * k2_Vy / 2) * sqrt(pow((Vx + h * k2_Vx / 2), 2) + pow((Vy + h * k2_Vy / 2), 2)) / 2;
        Ky_Cy3 = Ro(H) * Cya * alpha * S * (Vx + h * k2_Vx / 2) * sqrt(pow((Vx + h * k2_Vx / 2), 2) + pow((Vy + h * k2_Vy / 2), 2)) / 2;
        Ky_G3 = (Pi_0 * (ry + h * k2_ry / 2)) / pow(a3, 3);

        k3_rx = Vx + h * k2_Vx / 2;
        k3_ry = Vy + h * k2_Vy / 2;
        k3_Vx = (1 / m(t, m_0, m_k, t_k)) * (Kx_P3 - Kx_Cx3 - Kx_Cy3) - Kx_G3;
        k3_Vy = (1 / m(t, m_0, m_k, t_k)) * (Ky_P3 - Ky_Cx3 + Ky_Cy3) - Ky_G3;

        a4 = sqrt(pow(rx + h * k3_rx, 2) + pow(ry + h * k3_ry, 2));

        Kx_P4 = ((P_ud_p * m_k / t_k) * g) * cos(ny);
        Kx_Cx4 = Ro(H) * Cx0 * S * (Vx + h * k3_Vx) * sqrt(pow((Vx + h * k3_Vx), 2) + pow((Vy + h * k3_Vy), 2)) / 2;
        Kx_Cy4 = Ro(H) * Cya * alpha * S * (Vy + h * k3_Vy) * sqrt(pow((Vx + h * k3_Vx), 2) + pow((Vy + h * k3_Vy), 2)) / 2;
        Kx_G4 = (Pi_0 * (rx + h * k3_rx)) / pow(a4, 3);
        Ky_P4 = ((P_ud_p * m_k / t_k) * g) * sin(ny);
        Ky_Cx4 = Ro(H) * Cx0 * S * (Vy + h * k3_Vy) * sqrt(pow((Vx + h * k3_Vx), 2) + pow((Vy + h * k3_Vy), 2)) / 2;
        Ky_Cy4 = Ro(H) * Cya * alpha * S * (Vx + h * k3_Vx) * sqrt(pow((Vx + h * k3_Vx), 2) + pow((Vy + h * k3_Vy), 2)) / 2;
        Ky_G4 = (Pi_0 * (ry + h * k3_ry)) / pow(a4, 3);

        k4_rx = Vx + h * k3_Vx;
        k4_ry = Vy + h * k3_Vy;
        k4_Vx = (1 / m(t, m_0, m_k, t_k)) * (Kx_P4 - Kx_Cx4 - Kx_Cy4) - Kx_G4;
        k4_Vy = (1 / m(t, m_0, m_k, t_k)) * (Ky_P4 - Ky_Cx4 + Ky_Cy4) - Ky_G4;

        r_old[0] = rx, r_old[1] = ry;

        rx += (h * (k1_rx + 2 * k2_rx + 2 * k3_rx + k4_rx) / 6);
        ry += (h * (k1_ry + 2 * k2_ry + 2 * k3_ry + k4_ry) / 6);
        Vx += (h * (k1_Vx + 2 * k2_Vx + 2 * k3_Vx + k4_Vx) / 6);
        Vy += (h * (k1_Vy + 2 * k2_Vy + 2 * k3_Vy + k4_Vy) / 6);
        w1 += wd1 * h;

        del_r_old = del_r;

        r_new[0] = rx, r_new[1] = ry;

        del_rx = r_new[0] - r_new[0];

        del_ry = r_new[1] - r_old[1];

        del_r = sqrt(pow(del_rx, 2) + pow(del_ry, 2));

        S_old = del_r_old + del_r;

        r = sqrt(pow(rx, 2) + pow(ry, 2));

        V = sqrt(pow(Vx, 2) + pow(Vy, 2));

        H = (sqrt(pow(rx, 2) + pow(ry, 2)) - Rz);

        eta = atan(rx / ry);

        L = Rz * eta;

        alpha = ny - teta_v;

        q = Ro(H) * pow(V, 2) / 2;

        teta = ny + eta;

        P_not_cos_sin = ((P_ud_p * m_k / t_k) * g);

        peregruz_x = ((P_not_cos_sin + 2 * P_not_cos_sin + 2 * P_not_cos_sin + P_not_cos_sin) - (Kx_Cx1 + 2 * Kx_Cx2 + 2 * Kx_Cx3 + Kx_Cx4) - (Kx_Cy1 + 2 * Kx_Cy2 + 2 * Kx_Cy3 + Kx_Cy4)) / (m(t, m_0, m_k, t_k) * 9.8 * 6);

        peregruz_y = ((P_not_cos_sin + 2 * P_not_cos_sin + 2 * P_not_cos_sin + P_not_cos_sin) - (Ky_Cx1 + 2 * Ky_Cx2 + 2 * Ky_Cx3 + Ky_Cx4) + (Ky_Cy1 + 2 * Ky_Cy2 + 2 * Ky_Cy3 + Ky_Cy4)) / (m(t, m_0, m_k, t_k) * 9.8 * 6);
        }
        else
        {
        P = (P_ud_p * m_k / t_k) * g;
        Rx = 0;
        Ry = 0;
        wd1 = (P + sqrt(pow(Rx, 2) + pow(Ry, 2))) / m(t, m_0, m_k, t_k);

        a1 = sqrt(pow(rx, 2) + pow(ry, 2));
        Kx_P1 = ((P_ud_p * m_k / t_k) * g) * cos(ny);
        Kx_Cx1 = 0;
        Kx_Cy1 = 0;
        Kx_G1 = (Pi_0 * rx) / pow(a1, 3);
        Ky_P1 = ((P_ud_p * m_k / t_k) * g) * sin(ny);
        Ky_Cx1 = 0;
        Ky_Cy1 = 0;
        Ky_G1 = (Pi_0 * ry) / pow(a1, 3);

        k1_rx = Vx;
        k1_ry = Vy;
        k1_Vx = (1 / m(t, m_0, m_k, t_k)) * (Kx_P1 - Kx_Cx1 - Kx_Cy1) - Kx_G1;
        k1_Vy = (1 / m(t, m_0, m_k, t_k)) * (Ky_P1 - Ky_Cx1 + Ky_Cy1) - Ky_G1;

        a2 = sqrt(pow(rx + h * k1_rx / 2, 2) + pow(ry + h * k1_ry / 2, 2));
        Kx_P2 = ((P_ud_p * m_k / t_k) * g) * cos(ny);
        Kx_Cx2 = 0;
        Kx_Cy2 = 0;
        Kx_G2 = (Pi_0 * (rx + h * k1_rx / 2)) / pow(a2, 3);
        Ky_P2 = ((P_ud_p * m_k / t_k) * g) * sin(ny);
        Ky_Cx2 = 0;
        Ky_Cy2 = 0;
        Ky_G2 = (Pi_0 * (ry + h * k1_ry)) / pow(a2, 3);

        k2_rx = Vx + h * k1_Vx / 2;
        k2_ry = Vy + h * k1_Vy / 2;
        k2_Vx = (1 / m(t, m_0, m_k, t_k)) * (Kx_P2 - Kx_Cx2 - Kx_Cy2) - Kx_G2;
        k2_Vy = (1 / m(t, m_0, m_k, t_k)) * (Ky_P2 - Ky_Cx2 + Ky_Cy2) - Ky_G2;

        a3 = sqrt(pow(rx + h * k2_rx / 2, 2) + pow(ry + h * k2_ry / 2, 2));
        Kx_P3 = ((P_ud_p * m_k / t_k) * g) * cos(ny);
        Kx_Cx3 = 0;
        Kx_Cy3 = 0;
        Kx_G3 = (Pi_0 * (rx + h * k2_rx / 2)) / pow(a3, 3);
        Ky_P3 = ((P_ud_p * m_k / t_k) * g) * sin(ny);
        Ky_Cx3 = 0;
        Ky_Cy3 = 0;
        Ky_G3 = (Pi_0 * (ry + h * k2_ry / 2)) / pow(a3, 3);

        k3_rx = Vx + h * k2_Vx / 2;
        k3_ry = Vy + h * k2_Vy / 2;
        k3_Vx = (1 / m(t, m_0, m_k, t_k)) * (Kx_P3 - Kx_Cx3 - Kx_Cy3) - Kx_G3;
        k3_Vy = (1 / m(t, m_0, m_k, t_k)) * (Ky_P3 - Ky_Cx3 + Ky_Cy3) - Ky_G3;

        a4 = sqrt(pow(rx + h * k3_rx, 2) + pow(ry + h * k3_ry, 2));
        Kx_P4 = ((P_ud_p * m_k / t_k) * g) * cos(ny);
        Kx_Cx4 = 0;
        Kx_Cy4 = 0;
        Kx_G4 = (Pi_0 * (rx + h * k3_rx)) / pow(a4, 3);
        Ky_P4 = ((P_ud_p * m_k / t_k) * g) * sin(ny);
        Ky_Cx4 = 0;
        Ky_Cy4 = 0;
        Ky_G4 = (Pi_0 * (ry + h * k3_ry)) / pow(a4, 3);

        k4_rx = Vx + h * k3_Vx;
        k4_ry = Vy + h * k3_Vy;
        k4_Vx = (1 / m(t, m_0, m_k, t_k)) * (Kx_P4 - Kx_Cx4 - Kx_Cy4) - Kx_G4;
        k4_Vy = (1 / m(t, m_0, m_k, t_k)) * (Ky_P4 - Ky_Cx4 + Ky_Cy4) - Ky_G4;

        r_old[0] = rx, r_old[1] = ry;

        rx += (h * (k1_rx + 2 * k2_rx + 2 * k3_rx + k4_rx) / 6);
        ry += (h * (k1_ry + 2 * k2_ry + 2 * k3_ry + k4_ry) / 6);
        Vx += (h * (k1_Vx + 2 * k2_Vx + 2 * k3_Vx + k4_Vx) / 6);
        Vy += (h * (k1_Vy + 2 * k2_Vy + 2 * k3_Vy + k4_Vy) / 6);
        w1 += wd1 * h;

        del_r_old = del_r;

        r_new[0] = rx, r_new[1] = ry;

        del_rx = r_new[0] - r_new[0];

        del_ry = r_new[1] - r_old[1];

        del_r = sqrt(pow(del_rx, 2) + pow(del_ry, 2));

        S_old = del_r_old + del_r;

        r = sqrt(pow(rx, 2) + pow(ry, 2));

        V = sqrt(pow(Vx, 2) + pow(Vy, 2));

        H = (sqrt(pow(rx, 2) + pow(ry, 2)) - Rz);

        eta = atan(rx / ry);

        L = Rz * eta;

        alpha = ny - teta_v;

        q = Ro(H) * pow(V, 2) / 2;

        teta = ny + eta;

        P_not_cos_sin = ((P_ud_p * m_k / t_k) * g);

        peregruz_x = ((P_not_cos_sin + 2 * P_not_cos_sin + 2 * P_not_cos_sin + P_not_cos_sin) - (Kx_Cx1 + 2 * Kx_Cx2 + 2 * Kx_Cx3 + Kx_Cx4) - (Kx_Cy1 + 2 * Kx_Cy2 + 2 * Kx_Cy3 + Kx_Cy4)) / (m(t, m_0, m_k, t_k) * 9.8 * 6);

        peregruz_y = ((P_not_cos_sin + 2 * P_not_cos_sin + 2 * P_not_cos_sin + P_not_cos_sin) - (Ky_Cx1 + 2 * Ky_Cx2 + 2 * Ky_Cx3 + Ky_Cx4) + (Ky_Cy1 + 2 * Ky_Cy2 + 2 * Ky_Cy3 + Ky_Cy4)) / (m(t, m_0, m_k, t_k) * 9.8 * 6);
        }
        if (CheckBox1 -> Checked)
        NewSeries1 -> AddXY(t, H, "", OColor[ParamColor[1]]);
        if (CheckBox2 -> Checked)
        NewSeries2 -> AddXY(t, L, "", OColor[ParamColor[2]]);
        if (CheckBox3 -> Checked)
        NewSeries3 -> AddXY(t, alpha, "", OColor[ParamColor[3]]);
        if (CheckBox4 -> Checked)
        NewSeries4 -> AddXY(t, q, "", OColor[ParamColor[4]]);
        if (CheckBox5 -> Checked)
        NewSeries5 -> AddXY(t, V, "", OColor[ParamColor[5]]);
        if (CheckBox6 -> Checked)
        NewSeries6 -> AddXY(t, peregruz_x, "", OColor[ParamColor[6]]);
        if (CheckBox7 -> Checked)
        NewSeries7 -> AddXY(t, peregruz_y, "", OColor[ParamColor[7]]);
        if (CheckBox8 -> Checked)
        NewSeries8 -> AddXY(t, w1, "", OColor[ParamColor[8]]);
        if (CheckBox9 -> Checked)
        NewSeries9 -> AddXY(t, ny, "", OColor[ParamColor[9]]);
        if (CheckBox10 -> Checked)
        NewSeries10 -> AddXY(t, eta, "", OColor[ParamColor[10]]);
        if (CheckBox11 -> Checked)
        NewSeries11 -> AddXY(t, teta, "", OColor[ParamColor[11]]);
        output << t << setw(12) << w1 << setw(16) << V << setw(15) << q << setw(15)
               << peregruz_x << setw(14) << peregruz_y << setw(14) << ny << setw(15) << setw(14)
               << L << setw(15) << H << setw(15) << eta << setw(16) << teta << setw(16)
               << alpha << endl;
    }
    int m_const = m(t, m_0, m_k, t_k);
    t = 121;
    while(r > Rz)
    {
    t++;
        if (H > 100000)
        {
        ny = 0;
        P = 0;
        Rx = 0;
        Ry = 0;
        wd1 = (P + sqrt(pow(Rx, 2) + pow(Ry, 2))) / m_const;

        a1 = sqrt(pow(rx, 2) + pow(ry, 2));
        Kx_P1 = 0;
        Kx_Cx1 = 0;
        Kx_Cy1 = 0;
        Kx_G1 = (Pi_0 * rx) / pow(a1, 3);
        Ky_P1 = 0;
        Ky_Cx1 = 0;
        Ky_Cy1 = 0;
        Ky_G1 = (Pi_0 * ry) / pow(a1, 3);

        k1_rx = Vx;
        k1_ry = Vy;
        k1_Vx = (1 / m_const) * (Kx_P1 - Kx_Cx1 - Kx_Cy1) - Kx_G1;
        k1_Vy = (1 / m_const) * (Ky_P1 - Ky_Cx1 + Ky_Cy1) - Ky_G1;

        a2 = sqrt(pow(rx + h * k1_rx / 2, 2) + pow(ry + h * k1_ry / 2, 2));
        Kx_P2 = 0;
        Kx_Cx2 = 0;
        Kx_Cy2 = 0;
        Kx_G2 = (Pi_0 * (rx + h * k1_rx / 2)) / pow(a2, 3);
        Ky_P2 = 0;
        Ky_Cx2 = 0;
        Ky_Cy2 = 0;
        Ky_G2 = (Pi_0 * (ry + h * k1_ry)) / pow(a2, 3);

        k2_rx = Vx + h * k1_Vx / 2;
        k2_ry = Vy + h * k1_Vy / 2;
        k2_Vx = (1 / m_const) * (Kx_P2 - Kx_Cx2 - Kx_Cy2) - Kx_G2;
        k2_Vy = (1 / m_const) * (Ky_P2 - Ky_Cx2 + Ky_Cy2) - Ky_G2;

        a3 = sqrt(pow(rx + h * k2_rx / 2, 2) + pow(ry + h * k2_ry / 2, 2));
        Kx_P3 = 0;
        Kx_Cx3 = 0;
        Kx_Cy3 = 0;
        Kx_G3 = (Pi_0 * (rx + h * k2_rx / 2)) / pow(a3, 3);
        Ky_P3 = 0;
        Ky_Cx3 = 0;
        Ky_Cy3 = 0;
        Ky_G3 = (Pi_0 * (ry + h * k2_ry / 2)) / pow(a3, 3);

        k3_rx = Vx + h * k2_Vx / 2;
        k3_ry = Vy + h * k2_Vy / 2;
        k3_Vx = (1 / m_const) * (Kx_P3 - Kx_Cx3 - Kx_Cy3) - Kx_G3;
        k3_Vy = (1 / m_const) * (Ky_P3 - Ky_Cx3 + Ky_Cy3) - Ky_G3;

        a4 = sqrt(pow(rx + h * k3_rx, 2) + pow(ry + h * k3_ry, 2));
        Kx_P4 = 0;
        Kx_Cx4 = 0;
        Kx_Cy4 = 0;
        Kx_G4 = (Pi_0 * (rx + h * k3_rx)) / pow(a4, 3);
        Ky_P4 = 0;
        Ky_Cx4 = 0;
        Ky_Cy4 = 0;
        Ky_G4 = (Pi_0 * (ry + h * k3_ry)) / pow(a4, 3);

        k4_rx = Vx + h * k3_Vx;
        k4_ry = Vy + h * k3_Vy;
        k4_Vx = (1 / m_const) * (Kx_P4 - Kx_Cx4 - Kx_Cy4) - Kx_G4;
        k4_Vy = (1 / m_const) * (Ky_P4 - Ky_Cx4 + Ky_Cy4) - Ky_G4;

        r_old[0] = rx, r_old[1] = ry;

        rx += (h * (k1_rx + 2 * k2_rx + 2 * k3_rx + k4_rx) / 6);
        ry += (h * (k1_ry + 2 * k2_ry + 2 * k3_ry + k4_ry) / 6);
        Vx += (h * (k1_Vx + 2 * k2_Vx + 2 * k3_Vx + k4_Vx) / 6);
        Vy += (h * (k1_Vy + 2 * k2_Vy + 2 * k3_Vy + k4_Vy) / 6);
        w1 += wd1 * h;

        del_r_old = del_r;

        r_new[0] = rx, r_new[1] = ry;

        del_rx = r_new[0] - r_new[0];

        del_ry = r_new[1] - r_old[1];

        del_r = sqrt(pow(del_rx, 2) + pow(del_ry, 2));

        S_old = del_r_old + del_r;

        P_not_cos_sin = 0;

        teta_v = atan(Vy / Vx);

        V = sqrt(pow(Vx, 2) + pow(Vy, 2));

        H = (sqrt(pow(rx, 2) + pow(ry, 2)) - Rz);

        r = sqrt(pow(rx, 2) + pow(ry, 2));

        eta = atan(rx / ry);

        L = Rz * eta;

        alpha = 0;

        teta += eta;

        peregruz_x = ((P_not_cos_sin + 2 * P_not_cos_sin + 2 * P_not_cos_sin + P_not_cos_sin) - (Kx_Cx1 + 2 * Kx_Cx2 + 2 * Kx_Cx3 + Kx_Cx4) - (Kx_Cy1 + 2 * Kx_Cy2 + 2 * Kx_Cy3 + Kx_Cy4)) / (m_const * 9.8 * 6);

        peregruz_y = ((P_not_cos_sin + 2 * P_not_cos_sin + 2 * P_not_cos_sin + P_not_cos_sin) - (Ky_Cx1 + 2 * Ky_Cx2 + 2 * Ky_Cx3 + Ky_Cx4) + (Ky_Cy1 + 2 * Ky_Cy2 + 2 * Ky_Cy3 + Ky_Cy4)) / (m_const * 9.8 * 6);
        }
        else
        {
        ny = 0;
        P = 0;
        Rx = Cx0 * S * q;
        Ry = Cya * alpha * S * q;
        wd1 = (P + sqrt(pow(Rx, 2) + pow(Ry, 2))) / m_const;

        a1 = sqrt(pow(rx, 2) + pow(ry, 2));
        Kx_P1 = 0;
        Kx_Cx1 = Ro(H) * Cx0 * S * Vx * sqrt(pow(Vx, 2) + pow(Vy, 2)) / 2;
        Kx_Cy1 = Ro(H) * Cya * alpha * S * Vy * sqrt(pow(Vx, 2) + pow(Vy, 2)) / 2;
        Kx_G1 = (Pi_0 * rx) / pow(a1, 3);
        Ky_P1 = 0;
        Ky_Cx1 = Ro(H) * Cx0 * S * Vy * sqrt(pow(Vx, 2) + pow(Vy, 2)) / 2;
        Ky_Cy1 = Ro(H) * Cya * alpha * S * Vx * sqrt(pow(Vx, 2) + pow(Vy, 2)) / 2;
        Ky_G1 = (Pi_0 * ry) / pow(a1, 3);

        k1_rx = Vx;
        k1_ry = Vy;
        k1_Vx = (1 / m_const) * (Kx_P1 - Kx_Cx1 - Kx_Cy1) - Kx_G1;
        k1_Vy = (1 / m_const) * (Ky_P1 - Ky_Cx1 + Ky_Cy1) - Ky_G1;

        a2 = sqrt(pow(rx + h * k1_rx / 2, 2) + pow(ry + h * k1_ry / 2, 2));
        Kx_P2 = 0;
        Kx_Cx2 = Ro(H) * Cx0 * S * (Vx + h * k1_Vx / 2) * sqrt(pow((Vx + h * k1_Vx / 2), 2) + pow((Vy + h * k1_Vy / 2), 2)) / 2;
        Kx_Cy2 = Ro(H) * Cya * alpha * S * (Vy + h * k1_Vy / 2) * sqrt(pow((Vx + h * k1_Vx / 2), 2) + pow((Vy + h * k1_Vy / 2), 2)) / 2;
        Kx_G2 = (Pi_0 * (rx + h * k1_rx / 2)) / pow(a2, 3);
        Ky_P2 = 0;
        Ky_Cx2 = Ro(H) * Cx0 * S * (Vy + h * k1_Vy / 2) * sqrt(pow((Vx + h * k1_Vx / 2), 2) + pow((Vy + h * k1_Vy / 2), 2)) / 2;
        Ky_Cy2 = Ro(H) * Cya * alpha * S * (Vx + h * k1_Vx / 2) * sqrt(pow((Vx + h * k1_Vx / 2), 2) + pow((Vy + h * k1_Vy / 2), 2)) / 2;
        Ky_G2 = (Pi_0 * (ry + h * k1_ry)) / pow(a2, 3);

        k2_rx = Vx + h * k1_Vx / 2;
        k2_ry = Vy + h * k1_Vy / 2;
        k2_Vx = (1 / m_const) * (Kx_P2 - Kx_Cx2 - Kx_Cy2) - Kx_G2;
        k2_Vy = (1 / m_const) * (Ky_P2 - Ky_Cx2 + Ky_Cy2) - Ky_G2;

        a3 = sqrt(pow(rx + h * k2_rx / 2, 2) + pow(ry + h * k2_ry / 2, 2));
        Kx_P3 = 0;
        Kx_Cx3 = Ro(H) * Cx0 * S * (Vx + h * k2_Vx / 2) * sqrt(pow((Vx + h * k2_Vx / 2), 2) + pow((Vy + h * k2_Vy / 2), 2)) / 2;
        Kx_Cy3 = Ro(H) * Cya * alpha * S * (Vy + h * k2_Vy / 2) * sqrt(pow((Vx + h * k2_Vx / 2), 2) + pow((Vy + h * k2_Vy / 2), 2)) / 2;
        Kx_G3 = (Pi_0 * (rx + h * k2_rx / 2)) / pow(a3, 3);
        Ky_P3 = 0;
        Ky_Cx3 = Ro(H) * Cx0 * S * (Vy + h * k2_Vy / 2) * sqrt(pow((Vx + h * k2_Vx / 2), 2) + pow((Vy + h * k2_Vy / 2), 2)) / 2;
        Ky_Cy3 = Ro(H) * Cya * alpha * S * (Vx + h * k2_Vx / 2) * sqrt(pow((Vx + h * k2_Vx / 2), 2) + pow((Vy + h * k2_Vy / 2), 2)) / 2;
        Ky_G3 = (Pi_0 * (ry + h * k2_ry / 2)) / pow(a3, 3);

        k3_rx = Vx + h * k2_Vx / 2;
        k3_ry = Vy + h * k2_Vy / 2;
        k3_Vx = (1 / m_const) * (Kx_P3 - Kx_Cx3 - Kx_Cy3) - Kx_G3;
        k3_Vy = (1 / m_const) * (Ky_P3 - Ky_Cx3 + Ky_Cy3) - Ky_G3;

        a4 = sqrt(pow(rx + h * k3_rx, 2) + pow(ry + h * k3_ry, 2));
        Kx_P4 = 0;
        Kx_Cx4 = Ro(H) * Cx0 * S * (Vx + h * k3_Vx) * sqrt(pow((Vx + h * k3_Vx), 2) + pow((Vy + h * k3_Vy), 2)) / 2;
        Kx_Cy4 = Ro(H) * Cya * alpha * S * (Vy + h * k3_Vy) * sqrt(pow((Vx + h * k3_Vx), 2) + pow((Vy + h * k3_Vy), 2)) / 2;
        Kx_G4 = (Pi_0 * (rx + h * k3_rx)) / pow(a4, 3);
        Ky_P4 = 0;
        Ky_Cx4 = Ro(H) * Cx0 * S * (Vy + h * k3_Vy) * sqrt(pow((Vx + h * k3_Vx), 2) + pow((Vy + h * k3_Vy), 2)) / 2;
        Ky_Cy4 = Ro(H) * Cya * alpha * S * (Vx + h * k3_Vx) * sqrt(pow((Vx + h * k3_Vx), 2) + pow((Vy + h * k3_Vy), 2)) / 2;
        Ky_G4 = (Pi_0 * (ry + h * k3_ry)) / pow(a4, 3);

        k4_rx = Vx + h * k3_Vx;
        k4_ry = Vy + h * k3_Vy;
        k4_Vx = (1 / m_const) * (Kx_P4 - Kx_Cx4 - Kx_Cy4) - Kx_G4;
        k4_Vy = (1 / m_const) * (Ky_P4 - Ky_Cx4 + Ky_Cy4) - Ky_G4;

        r_old[0] = rx, r_old[1] = ry;

        rx += (h * (k1_rx + 2 * k2_rx + 2 * k3_rx + k4_rx) / 6);
        ry += (h * (k1_ry + 2 * k2_ry + 2 * k3_ry + k4_ry) / 6);
        Vx += (h * (k1_Vx + 2 * k2_Vx + 2 * k3_Vx + k4_Vx) / 6);
        Vy += (h * (k1_Vy + 2 * k2_Vy + 2 * k3_Vy + k4_Vy) / 6);
        w1 += wd1 * h;

        del_r_old = del_r;

        r_new[0] = rx, r_new[1] = ry;

        del_rx = r_new[0] - r_new[0];

        del_ry = r_new[1] - r_old[1];

        del_r = sqrt(pow(del_rx, 2) + pow(del_ry, 2));

        S_old = del_r_old + del_r;

        teta_v = atan(Vy / Vx);

        V = sqrt(pow(Vx, 2) + pow(Vy, 2));

        H = (sqrt(pow(rx, 2) + pow(ry, 2)) - Rz);

        r = sqrt(pow(rx, 2) + pow(ry, 2));

        eta = atan(rx / ry);

        L = Rz * eta;

        alpha = 0;

        teta += eta;

        peregruz_x = ((P_not_cos_sin + 2 * P_not_cos_sin + 2 * P_not_cos_sin + P_not_cos_sin) - (Kx_Cx1 + 2 * Kx_Cx2 + 2 * Kx_Cx3 + Kx_Cx4) - (Kx_Cy1 + 2 * Kx_Cy2 + 2 * Kx_Cy3 + Kx_Cy4)) / (m_const * 9.8 * 6);

        peregruz_y = ((P_not_cos_sin + 2 * P_not_cos_sin + 2 * P_not_cos_sin + P_not_cos_sin) - (Ky_Cx1 + 2 * Ky_Cx2 + 2 * Ky_Cx3 + Ky_Cx4) + (Ky_Cy1 + 2 * Ky_Cy2 + 2 * Ky_Cy3 + Ky_Cy4)) / (m_const * 9.8 * 6);
        }
        if (CheckBox1 -> Checked)
        NewSeries1 -> AddXY(t, H, "", OColor[ParamColor[1]]);
        if (CheckBox2 -> Checked)
        NewSeries2 -> AddXY(t, L, "", OColor[ParamColor[2]]);
        if (CheckBox3 -> Checked)
        NewSeries3 -> AddXY(t, alpha, "", OColor[ParamColor[3]]);
        if (CheckBox4 -> Checked)
        NewSeries4 -> AddXY(t, q, "", OColor[ParamColor[4]]);
        if (CheckBox5 -> Checked)
        NewSeries5 -> AddXY(t, V, "", OColor[ParamColor[5]]);
        if (CheckBox6 -> Checked)
        NewSeries6 -> AddXY(t, peregruz_x, "", OColor[ParamColor[6]]);
        if (CheckBox7 -> Checked)
        NewSeries7 -> AddXY(t, peregruz_y, "", OColor[ParamColor[7]]);
        if (CheckBox8 -> Checked)
        NewSeries8 -> AddXY(t, w1, "", OColor[ParamColor[8]]);
        if (CheckBox9 -> Checked)
        NewSeries9 -> AddXY(t, ny, "", OColor[ParamColor[9]]);
        if (CheckBox10 -> Checked)
        NewSeries10 -> AddXY(t, eta, "", OColor[ParamColor[10]]);
        if (CheckBox11 -> Checked)
        NewSeries11 -> AddXY(t, teta, "", OColor[ParamColor[11]]);
        output << t << setw(12) << w1 << setw(16) << V << setw(15) << q << setw(15)
               << peregruz_x << setw(14) << peregruz_y << setw(14) << ny << setw(15) << setw(14)
               << L << setw(15) << H << setw(15) << eta << setw(16) << teta << setw(16)
               << alpha << endl;
    }
}
//---------------------------------------------------------------------------
void __fastcall TForm1::Button2Click(TObject *Sender)
{
   Chart1 -> RemoveAllSeries();
   Label11 -> Visible = false;
   Label12 -> Visible = false;
   Label13 -> Visible = false;
}
//---------------------------------------------------------------------------
void __fastcall TForm1::Chart1MouseDown(TObject *Sender,
      TMouseButton Button, TShiftState Shift, int X, int Y)
{
     Label13 -> Left = X;
     Label13 -> Top = Y;
     if (CheckBox1 -> Checked)
     {
       Screen -> Cursor = crCross;
       NewSeries1 -> GetCursorValues(t, H);
       Label13 -> Visible = true;
       Label13 -> Caption = "      H = " + String(floor(H)) + " м" + "    t = " + String(floor(t)) + " с";
     }
     else
     if (CheckBox2 -> Checked)
     {
       Screen -> Cursor = crCross;
       NewSeries2 -> GetCursorValues(t, L);
       Label13 -> Visible = true;
       Label13 -> Caption = "      L = " + String(floor(L)) + " м" + "    t = " + String(floor(t)) + " с";
     }
}
//---------------------------------------------------------------------------
void __fastcall TForm1::PageControl1MouseMove(TObject *Sender,
      TShiftState Shift, int X, int Y)
{
   Screen -> Cursor = crDefault;  
}
//---------------------------------------------------------------------------
