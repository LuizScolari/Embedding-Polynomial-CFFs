#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/ZZ_pXFactoring.h>
#include <vector>
#include <iostream>
#include <unordered_map>     
#include <optional>   
#include <utility>
#include <sstream>
#include <functional>

using namespace NTL;

struct ZZ_pE_Hash {
    size_t operator()(const ZZ_pE& x) const {
        std::ostringstream oss;
        oss << x;
        return std::hash<std::string>{}(oss.str());
    }
};

class ExtensionElementsGenerator {
private:
    std::vector<ZZ_pE> elems;
    long p;
    int n;
    
public:
    ExtensionElementsGenerator(long prime, int degree)
        : p(prime), n(degree)
    {
        // 1) Inicializa o corpo base F_p
        ZZ_p::init(conv<ZZ>(p));
        
        // 2) Gera um polinômio irreduzível P(x) de grau n em F_p[x]
        ZZ_pX P;
        BuildIrred(P, n);
        
        // 3) Inicializa o corpo de extensão F_{p^n}
        ZZ_pE::init(P);
        
        // 4) Calcula q = p^n e reserva espaço
        long q = 1;
        for (int i = 0; i < n; ++i) q *= p;
        elems.reserve(q);
        
        // 5) Enumera cada idx de 0 a p^n−1
        for (long idx = 0; idx < q; ++idx) {
            long tmp = idx;
            ZZ_pX f; // polinômio representando o elemento
            for (int d = 0; d < n; ++d) {
                long coef = tmp % p;
                if (coef != 0)
                    SetCoeff(f, d, coef);
                tmp /= p;
            }
            ZZ_pE e; 
            conv(e, f);
            elems.push_back(e);
        }
    }
    
    const std::vector<ZZ_pE>& getElements() const {
        return elems;
    }
    
    void printElements() const {
        std::cout << "Extension field GF(" << p << "^" << n << ") elements:\n";
        long q = elems.size();
        for (long i = 0; i < q; ++i) {
            ZZ_pX poly = rep(elems[i]); // rep() traz o polinômio interno
            // imprime sempre n coeficientes, de x^0 a x^(n-1)
            std::cout << "Element " << i << ": [";
            for (int d = 0; d < n; ++d) {
                // coeff(poly, d) é um ZZ_p, imprime como inteiro
                std::cout << coeff(poly, d);
                if (d+1 < n) std::cout << " ";
            }
            std::cout << "]\n";
        }
    }
};

std::vector<ZZ_pEX>
generate_polynomials(const std::vector<ZZ_pE>& elems, int k) {
    long q = elems.size();
    // total de polinômios de grau ≤ k = q^(k+1)
    long total = 1;
    for(int i = 0; i <= k; ++i)
        total *= q;
    
    std::vector<ZZ_pEX> polys;
    polys.reserve(total);
    
    // agora varremos de 0 até total-1
    for(long idx = 0; idx < total; ++idx){
        long tmp = idx;
        ZZ_pEX f;
        for(int i = 0; i <= k; ++i){
            long ci = tmp % q;
            SetCoeff(f, i, elems[ci]);
            tmp /= q;
        }
        polys.push_back(f);
    }
    return polys;
}

std::vector<std::pair<ZZ_pE, ZZ_pE>> 
generate_rows(const std::vector<ZZ_pE>& elems){
    std::vector<std::pair<ZZ_pE, ZZ_pE>> tuplas;
    tuplas.reserve(elems.size() * elems.size());
    for(std::size_t i=0; i<elems.size(); i++){
        for(std::size_t j=0; j<elems.size(); j++){
            tuplas.push_back({elems[i], elems[j]});
        }
    }
    return tuplas;
}

class HashMapVetores {
private:
    std::unordered_map<ZZ_pE, std::vector<ZZ_pE>, ZZ_pE_Hash> hash_map;

public:
    HashMapVetores(const std::vector<ZZ_pE>& chaves, const std::vector<std::vector<ZZ_pE>>& vetores) {
        adicionar(chaves, vetores);
    }
    
    HashMapVetores() = default;
    
    void adicionar(const std::vector<ZZ_pE>& chave, const std::vector<std::vector<ZZ_pE>>& vetorE) {
        for(size_t i = 0; i < chave.size(); i++) {
            hash_map[chave[i]] = vetorE[i];
        }
    }
    
    std::optional<std::vector<ZZ_pE>> obter(const ZZ_pE& chave) const {
        auto it = hash_map.find(chave);
        if (it != hash_map.end()) {
            return it->second;
        }
        return std::nullopt;
    }
};

void generate_cffs(int p, int n, int k){

    ExtensionElementsGenerator gen(p, n);
    std::vector<ZZ_pE> elems = gen.getElements();
    auto polys = generate_polynomials(elems, k);

    std::vector<std::vector<ZZ_pE>> results;
    results.reserve(elems.size());
    for (const auto& xi : elems) {
        std::vector<ZZ_pE> row;
        row.reserve(polys.size());
        for (const auto& P : polys) {
            // eval(polinômio, ponto)
            row.push_back(eval(P, xi));
        }
        results.push_back(std::move(row));
    }

    HashMapVetores hash_map(elems, results);
    auto rows = generate_rows(elems);

    std::vector<std::vector<int>> cff;
    for (const auto& [a, b] : rows) {
        std::vector<int> row;
        
        auto opt_row_E = hash_map.obter(a);
        if (opt_row_E.has_value()) {
            std::vector<ZZ_pE> row_E = opt_row_E.value();
            
            for(size_t i = 0; i < row_E.size(); i++) {
                if(b == row_E[i]) {
                    row.push_back(1);
                } else {
                    row.push_back(0);
                }
            }
        }
        cff.push_back(row);
    }

    for (size_t i = 0; i < cff.size(); i++) {
        for (int val : cff[i]) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }
}

/*
void test(){
    int p = 2;
    int n = 2;
    int k = 1;

    ExtensionElementsGenerator gen(p, n);
    auto polys = generate_polynomials(gen.getElements(), k);
    auto mat = evaluate_Polynomials(p, n, k);
    
    for (std::size_t i = 0; i < mat.size(); ++i) {
        for (std::size_t j = 0; j < mat[i].size(); ++j) {
            std::cout << mat[i][j] << " ";
        }
        std::cout << "\n";
    }

    auto rows = generate_rows(gen.getElements());

    for (const auto& [a, b] : rows) {
        std::cout << "(" << a << ", " << b << ")" << std::endl;
    } 
}
*/

int main() {
    int p = 2;
    int n = 2;
    int k = 1;

    generate_cffs(p, n, k);  
    
    return 0;
}