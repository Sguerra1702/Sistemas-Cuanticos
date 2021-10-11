import numpy as np


def probparticleinline():
    print("dimensión de su estado cuantico")
    dim = int(input())
    npket = [0 for rep in range(dim)]
    for rep in range(dim):
        print("ingrese el ", rep +1, "componente de su vector")
        numero = complex(input())
        npket[rep] = numero
    print(npket)
    ket = np.array(npket)
    print(ket)
    print("De cual posicion desea hallar la probabilidad")
    pos = int(input())
    while pos > dim:
        print("No es posible calcular esta posición, intente de nuevo")
        pos = int(input())
    print(ket[pos])
    print(abs(ket[pos]))
    prob_punto_linea = abs(ket[pos])/np.linalg.norm(ket)
    print("la probabilidad de encontrar al punto en la posición", pos, "es:", prob_punto_linea)
    return prob_punto_linea, dim, ket


def prog_drill_421():
    print("Cual es la dimensión de su observable y estado")
    dim = int(input())
    estado = [0 for rep in range(dim)]
    for n in range(dim):
        print("Digite la", n, "° posición de su estado")
        estado[n] = complex(input())
    ket = np.array(estado)
    obs = np.zeros((dim, dim), dtype = complex)
    conj = np.zeros((dim, dim), dtype = complex)
    for i in range(dim):
        for j in range(dim):
            print("digite el ", i, "° dato de su observable")
            obs[i, j] = complex(input())
            conj[i, j] = np.conjugate(obs[i, j])
    conj_tras = np.transpose(conj)
    print("el observable en el sistema es:")
    print(obs)
    print("La adjunta de la matriz es")
    print(conj_tras)
    check_hermitian = 0
    for i in range(dim):
        for j in range(dim):
            if obs[i, j] == conj_tras[i, j]:
                check_hermitian += 1
    if check_hermitian == dim**2:
        omega_psi = obs.dot(ket)
        print("\n omega psi =")
        print(omega_psi)
        med_value = float(np.dot(np.conjugate(omega_psi), ket))
        print("\nel valor medio del observable es:")
        print(med_value)
        miu = np.multiply(med_value, np.identity(dim))
        print("\n la matriz valor medio es")
        print(miu)
        delta_miu = obs - miu
        print("\n el observable menos la matriz valor medio es")
        print(delta_miu)
        prod = np.matmul(delta_miu, delta_miu)
        print("\n delta miu ^2 =",prod)
        varianza = float((np.conjugate(ket).dot(prod)).dot(ket))
        print("\nLa varianza en el sistema cuántico es: ")
        print(varianza)
    else:
        print("La matriz del observable no es hermitiana, no puedo realizar las operaciones")


#Ejercicio 4.4.1
def ejercicio_441():
    u1 = [[0,1],
           [1, 0]]
    u_1 = np.asmatrix(u1)
    mU1 = np.matmul(u_1, np.transpose(u_1))
    u2 = [[(np.sqrt(2))/2, (np.sqrt(2))/2],
          [(np.sqrt(2))/2, -(np.sqrt(2))/2]]
    u_2 = np.asmatrix(u2)
    mU2 = np.matmul(u_2, np.transpose(u_2))
    I2 = np.identity(2)
    cont_unit_1 = 0
    cont_unit_2 = 0
    for i in range(len(u_1)):
        for j in range(len(u_2)):
            if mU1[i, j] == I2[i, j]:
                cont_unit_1 += 1
            if mU2[i, j] == I2[i, j]:
                cont_unit_2 += 2
    if cont_unit_1 == len(u_1)**2:
        print("La matriz:")
        print(mU1)
        print("es unitaria")
    else:
        print("La matriz:")
        print(mU1)
        print("no es unitaria")
    if cont_unit_2 == len(u_2)**2:
        print("La matriz:")
        print(mU2)
        print("es unitaria")
    else:
        print("La matriz:")
        print(mU2)
        print("no es unitaria")
    prod = np.matmul(u_1, u_2)
    cont_prod = 0
    for i in range(len(prod)):
        for j in range(len(prod)):
            if prod[i, j] == I2[i, j]:
                cont_prod += 1
    if cont_prod == len(prod)**2:
        print("La matriz:")
        print(u_2)
        print("es unitaria")
    else:
        print("La matriz:")
        print(u_2)
        print("no es unitaria")


#Ejercicio 4.4.1
def ejercicio_442():
    estado = [1, 0, 0, 0]
    ket = np.array(estado)
    mat = [[0, 1/(np.sqrt(2)), 1/(np.sqrt(2)), 0],
           [(1j/(np.sqrt(2))), 0, 0, 1/(np.sqrt(2))],
           [1/(np.sqrt(2)), 0, 0, (1j/(np.sqrt(2)))],
           [0, 1/(np.sqrt(2)), -1/(np.sqrt(2)), 0]]
    mapa_unit = np.asmatrix(mat)
    pot = np.linalg.matrix_power(mapa_unit, 3)
    print(pot)
    fin_estado = pot.dot(ket)
    print("la condicion final del estado es:")
    print(fin_estado)
    print(ket[4])
    print(abs(ket[4]))
    prob_punto_linea = abs(ket[4])/np.linalg.norm(ket)
    print("la probabilidad de encontrar al punto en la posición 3 es:", prob_punto_linea)


def main():
    print("Librería de sistemas cuánticos capítulo 4")
    print("\nEjercicio Programming drill 4.1.1")
    probparticleinline()
    print("\nEjercicio Programming drill 4.2.1")
    prog_drill_421()


main()
