#ifndef PHYSICCONST_H
#define PHYSICCONST_H
#define PI 3.14159265358979323846264338327950288419716939937510 
#define CC 2.99792458e+8 //�������� ����� � �������
#define CCC 2.99792458e+10 //�������� ����� � �������
#define MU_Z 0.0000012566370614359173//���������� ��������� ��-�� � �������  4.0*M_PI*1.0e-7 
#define EPS_Z 0.0000000000088541878176203892 //��������������� ������������� 1.0/(CC*CC*MU_Z)
#define FREQ 5.0e+9 // ������� ���������
#define LAMBDA (CC / FREQ) // ���� ������ ���� �� ������ ��� ���� ������� ����� ������ �������
#define OMEGA (2.0 * PI * FREQ) // �� �� �������
#define IMP_Z sqrt(MU_Z / EPS_Z)

#endif