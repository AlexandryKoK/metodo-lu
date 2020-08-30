program metodolu

!Universidade de Brasilia - Metodos Computacionais A
!Professor: Luiz A. Ribeiro Junior - Instituto de Física
!Aluno: Alexandry Moreira Alves Pinto - 17/0078761
!Problema 1 - Solução de sistema de equaçoẽs lineares usando metodo LU.
!Github url: 

integer*8 :: n, k
real*8 :: temoum, termodois,termoy, termox
double precision, dimension(100,100) :: matrix
double precision, dimension(100,100) :: matrixL
double precision, dimension(100,100) :: matrixU
double precision, dimension(100,1) :: matrixb
double precision, dimension(100,1) :: matrixy
double precision, dimension(100,1) :: matrixx



!O algoritmo é: temos Ax=b, decompomos a matriz A em uma operação de duas, LU
!reescrevemos como LUx=b, criamos o seguinte sistema, Ux=y -> Ly=b
!resolvendo para y, podemos achar x.


write(*,*) "Bem vindo ao programa de solucao de sistemas lineares usando metodo LU."
write(*,*) 
write(*,*) "Digite a ordem da sua matriz (numero inteiro e sem contar a coluna das variaveis):"
write(*,*)
read (*,*) n
write(*,*) "Quando for escrever fração, transforme em numero decimal."
write(*,*)


do i = 1, n
	do j = 1, n+1
		if (j < (n+1)) then
		write(*,*) "Digite o",i, "elemento da",j, "coluna"
		read(*,*) matrix(i,j)
		else 
		write(*,*) "Elemento",i, " da coluna de resultados."
		read(*,*) matrix(i,j)
		end if
	end do
end do

!Criar uma matriz b coluna onde os elementos são resultados das equações
do i=1,n
matrixb(i,1) = matrix(i,n+1)
end do

!1 Passo - Calculo da primeira coluna de [matrixL]
!matrixL(i,1) = matrix(i,1)

do i = 1,n

matrixL(i,1) = matrix(i,1)

end do



!2 Passo - Toda diagonal de U ser 1

do i=1,n
do j=1,n

if (i == j) then
matrixU(i,j) = 1
end if

end do
end do



!3 Passo - Calculo de elementos da primeira linha de U

do j=2,n

matrixU(1,j) = matrix(1,j) / matrixL(1,1)

end do




!4 Passo - Calculo dos demais elementos de L,U

do i=2,n

!Calculo dos elementos de L
do j=2,i
termoum = 0 !sempre que j for mudar temos que zerar a somatoria termoum


do k=1,j-1
termoum = termoum + matrixL(i,k) * matrixU(k,j)
end do


matrixL(i,j) = (matrix(i,j) - termoum) 
end do

!Calculo dos elementos de U
do j=i+1,n
termodois = 0

do k=1,i-1
termodois = termodois + matrixL(i,k) * matrixU(k,j)
end do

matrixU(i,j) = (matrix(i,j) - termodois) /  matrixL(i,i)

end do
end do


!5 Passo - [matrixL]*[matrixy] = [matrixb]; precisamos agora resolver para y

do i=1,n

if (i == 1) then

matrixy(i,1) = matrixb(i,1) / matrixL(i,i)
!write(*,*)matrixy(i,1)

else if (i /= 1) then
termoy = 0

do j=1,i-1
termoy = termoy + matrixL(i,j) * matrixy(j,1)
end do

matrixy(i,1) = (matrixb(i,1) - termoy) / matrixL(i,i)
!write(*,*) matrixy(i,1)
end if


end do


!6 Passo - [matrixU]*[matrixx] = [matrixy]; precisamos resolver para x

do i=n,1,-1

if (i==n) then
matrixx(i,1) = matrixy(i,1) / 1
else if (i /= n) then

termox = 0
do j=i+1,n
termox = termox +  matrixU(i,j) * matrixx(j,1)
end do

matrixx(i,1) = (matrixy(i,1) - termox)
end if

end do


!Final - Imprimir o valor das variaveis finais.

do i=1,n

write(*,*) "Variavel",i,"=", matrixx(i,1)

end do



!Se quiser imprimir os elementos de L e U
!do i=1,n
!do j=1,n
!write(*,*) matrixL(i,j)
!end do
!end do
!do i=1,n
!do j=1,n
!write(*,*) matrixU(i,j)
!end do
!end do

end program Metodolu