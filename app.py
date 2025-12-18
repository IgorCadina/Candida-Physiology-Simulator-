
##===================================================================#
##                                                                   #
##                 Candida Physiology Simulator                      #
##                                                                   #
##===================================================================#

import os
import base64
from io import BytesIO
import numpy as np
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['axes.unicode_minus'] = False
import matplotlib.pyplot as plt
from flask import Flask, request, render_template

app = Flask(__name__)

class SimuladorCandida:
    def __init__(self):
        self.meios = {
            'Saboroud Dextrose Líquido (SDB)': {
                'fator_crescimento': 1.1,
                'fator_hifa': 0.4,      
                'od_max': 14.0,
                'descricao': "Meio padrão para cultivo de fungos",
                'ph_otimo': 6.5,
                'cor': '#3B82F6'
            },
            'Saboroud Dextrose Ágar (SDA)': {
                'fator_crescimento': 1.1,
                'fator_hifa': 0.4,  
                'od_max': 14.0,
                'descricao': "Meio sólido com ágar 2%",
                'ph_otimo': 6.5,
                'cor': '#10B981'
            },
            'YPD (Yeast Extract Peptone Dextrose)': {
                'fator_crescimento': 1.3,
                'fator_hifa': 0.3, 
                'od_max': 18.0,
                'descricao': "Meio rico para crescimento de levedo",
                'ph_otimo': 6.0,
                'cor': '#F59E0B'
            },
            'Meio de Lee': {
                'fator_crescimento': 0.8,
                'fator_hifa': 1.0,  
                'od_max': 10.0,
                'descricao': "Meio otimizado para indução de filamentação",
                'ph_otimo': 7.0,
                'cor': '#EF4444'
            },
            'RPMI-1640': {
                'fator_crescimento': 0.7,
                'fator_hifa': 0.7,  
                'od_max': 8.0,
                'descricao': "Simula condições do hospedeiro humano",
                'ph_otimo': 7.4,
                'cor': '#8B5CF6'
            },
            'Meio Sintético Mínimo (SD)': {
                'fator_crescimento': 0.5,
                'fator_hifa': 0.2,  
                'od_max': 6.0,
                'descricao': "Meio mínimo definido",
                'ph_otimo': 6.0,
                'cor': '#06B6D4'
            },
        }
        
        self.parametros = {
            'ph': {
                'otimo': 6.5,
                'min': 2.0,
                'max': 9.0,
                'letal_min': 1.0,
                'letal_max': 12.0,
                'sigma': 1.5
            },
            'temperatura': {
                'otimo': 37.0,
                'min': 4.0,
                'max': 45.0,
                'letal': 55.0,
                'sigma': 5.0  
            },
            'farnesol': {
                'inibicao_hifas': 10.0,
                'inibicao_crescimento': 50.0,
                'limiar_toxicidade': 100.0
            },
            'crescimento': {
                'taxa_max': 0.46,
                'carrying_capacity_factor': 0.8
            }
        }

    def resposta_ph_hifas(self, ph):
        if ph >= 7.5:
            return 1.0
        elif ph >= 7.0:
            return 0.7 + 0.3 * ((ph - 7.0) / 0.5)
        elif ph >= 6.5:
            return 0.3 + 0.4 * ((ph - 6.5) / 0.5)
        elif ph >= 6.0:
            return 0.1 + 0.2 * ((ph - 6.0) / 0.5)
        else:
            return 0.05 * (ph / 6.0)
    
    def resposta_temp_hifas(self, temp):
        if temp >= 37.0:
            return 1.0
        elif temp >= 35.0:
            return 0.8 + 0.2 * ((temp - 35.0) / 2.0)
        elif temp >= 33.0:
            return 0.6 + 0.2 * ((temp - 33.0) / 2.0)
        elif temp >= 30.0:
            return 0.3 + 0.3 * ((temp - 30.0) / 3.0)
        elif temp >= 25.0:
            return 0.1 + 0.2 * ((temp - 25.0) / 5.0)
        else:
            return 0.05 * (temp / 25.0)

    def prever(self, meio, ph, temp, farnesol):
        meio_info = self.meios.get(meio, self.meios['Saboroud Dextrose Líquido (SDB)'])
        
        ph_otimo_meio = meio_info['ph_otimo']
        eficiencia_ph = np.exp(-((ph - ph_otimo_meio) ** 2) / (2 * self.parametros['ph']['sigma'] ** 2))
    
        eficiencia_temp = np.exp(-((temp - self.parametros['temperatura']['otimo']) ** 2) / 
                                 (2 * self.parametros['temperatura']['sigma'] ** 2))
        
        taxa_crescimento_basal = meio_info['fator_crescimento'] * eficiencia_temp * eficiencia_ph
        
        inibicao_crescimento = 0.0
        if farnesol > self.parametros['farnesol']['limiar_toxicidade']:
            inibicao_crescimento = min(0.7, (farnesol - 100.0) / 200.0)
        elif farnesol > self.parametros['farnesol']['inibicao_crescimento']:
            inibicao_crescimento = min(0.3, (farnesol - 50.0) / 100.0)
        
        taxa_crescimento_basal *= (1 - inibicao_crescimento)

        if ph >= 4.0 and ph <= 8.5 and temp >= 20.0 and temp <= 42.0:
            taxa_crescimento_basal = max(taxa_crescimento_basal, 0.3)
        
        glicolise = taxa_crescimento_basal
        ciclo_tca = taxa_crescimento_basal * 0.7
        atividade_ribossomal = taxa_crescimento_basal * 0.9
        
        induz_ph = self.resposta_ph_hifas(ph)
        induz_temp = self.resposta_temp_hifas(temp)
        
        formacao_hifas_base = meio_info['fator_hifa'] * induz_ph * induz_temp
        
        if induz_ph > 0.7 and induz_temp > 0.7:
            fator_sinergia = 1.0 + (induz_ph * induz_temp) * 0.3
            formacao_hifas_base *= fator_sinergia
        
        if meio == 'Meio de Lee' and ph >= 7.0 and temp >= 35.0:
            formacao_hifas_base = min(1.0, formacao_hifas_base * 1.2)
        
        if meio == 'RPMI-1640' and ph >= 7.2 and temp >= 36.0:
            formacao_hifas_base = min(1.0, formacao_hifas_base * 1.15)
        
        if farnesol > 0:
            K = self.parametros['farnesol']['inibicao_hifas'] 
            n = 2
            inibicao = (farnesol ** n) / (K ** n + farnesol ** n)
            formacao_hifas_base *= (1 - inibicao * 0.9)
        
        formacao_hifas_base = np.clip(formacao_hifas_base, 0.0, 1.0)
        
        estresse_total = 0.0
        
        if temp > 42.0:
            estresse_total += 0.7
        elif temp > 39.0:
            estresse_total += 0.4
        elif temp > 37.5:
            estresse_total += 0.1
        elif temp < 25.0:
            estresse_total += 0.3
        
        if ph < 4.0 or ph > 8.5:
            estresse_total += 0.6
        elif ph < 5.0 or ph > 8.0:
            estresse_total += 0.3
        
        ativacao_hog1 = np.clip(estresse_total, 0.0, 0.95)
        
        if ativacao_hog1 > 0.5:
            atividade_ribossomal *= 0.7
            parede_celular = 0.9
            resistencia_base = ativacao_hog1 * 0.8
        elif ativacao_hog1 > 0.3:
            atividade_ribossomal *= 0.85
            parede_celular = 0.75
            resistencia_base = ativacao_hog1 * 0.6
        else:
            parede_celular = 0.5 + (atividade_ribossomal * 0.3)
            resistencia_base = 0.1
        
        scores = {
            'Glicólise': np.clip(glicolise, 0, 1),
            'Ciclo TCA': np.clip(ciclo_tca, 0, 1),
            'Atividade Ribossomal': np.clip(atividade_ribossomal, 0, 1),
            'Formação de Hifas': formacao_hifas_base,
            'Biossíntese de Ergosterol': np.clip(atividade_ribossomal * 0.85, 0, 1),
            'Parede Celular': np.clip(parede_celular, 0, 1),
            'Via MAPK Hog1': ativacao_hog1,
            'Resistência a Azóis': np.clip(resistencia_base + parede_celular * 0.2, 0.1, 1),
            'Ácidos Graxos': np.clip(ciclo_tca * 0.8, 0, 1),
            'Consumo de Glicose': np.clip(glicolise * 0.9, 0, 1)
        }
        
        return scores, meio_info['od_max']


simulador = SimuladorCandida()

simulador = SimuladorCandida()

def gerar_grafico(ribo, od_inicial, horas, od_max):
    od_capacidade = od_max * (0.2 + 0.8 * ribo)
    
    if ribo > 0.01:
        taxa = simulador.parametros['crescimento']['taxa_max'] * ribo
    else:
        taxa = 0.001
    
    t = np.linspace(0, horas, 100)
    
    if ribo < 0.05:
        od = od_inicial * np.exp(-0.01 * t)
        estado = "Latência/Declínio"
        cor = '#94A3B8'
    else:
        a = (od_capacidade - od_inicial) / od_inicial
        od = od_capacidade / (1 + a * np.exp(-taxa * t))
        
        if ribo > 0.7:
            estado = "Crescimento Exponencial"
            cor = '#10B981'
        elif ribo > 0.3:
            estado = "Crescimento Linear"
            cor = '#F59E0B'
        else:
            estado = "Crescimento Lento"
            cor = '#EAB308'
    
    plt.style.use('seaborn-v0_8-whitegrid')
    fig, ax = plt.subplots(figsize=(12, 7), facecolor='#F9FAFB')
    
    ax.plot(t, od, color=cor, linewidth=3.5, label=f'{estado}')
    ax.fill_between(t, 0, od, alpha=0.15, color=cor)
    
    ax.axhline(y=od_capacidade, color='#EF4444', linestyle='--', 
               linewidth=2, alpha=0.7, label=f'Capacidade: {od_capacidade:.2f} OD')
    ax.axhline(y=od_inicial, color='#6B7280', linestyle=':', 
               linewidth=1.5, alpha=0.5, label=f'Inicial: {od_inicial:.2f} OD')
    
    if ribo > 0.3 and taxa > 0.1:
        t_inflexao = np.log(a) / taxa
        if 0 < t_inflexao < horas:
            ax.axvline(x=t_inflexao, color='#8B5CF6', linestyle='--', 
                       alpha=0.5, linewidth=1.5)
            ax.text(t_inflexao + 0.5, od_capacidade * 0.1, 
                    'Ponto de\ninflexão', fontsize=9, color='#8B5CF6',
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    ax.set_title(f"Cinética de Crescimento - {int(horas)} horas", 
                 fontsize=18, fontweight='bold', pad=20, color='#1F2937')
    ax.set_xlabel("Tempo (horas)", fontsize=14, fontweight='medium', color='#374151')
    ax.set_ylabel("OD600", fontsize=14, fontweight='medium', color='#374151')
    
    ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
    ax.set_xlim(0, horas)
    ax.set_ylim(0, od_capacidade * 1.15)
    
    ax.legend(loc='upper left', fontsize=11, frameon=True, 
              fancybox=True, shadow=True, facecolor='white')
    
    texto_info = f"""μ = {taxa:.3f} h-1
t_d = {0.693/taxa:.1f} h (se μ>0.01)
OD_final = {od[-1]:.3f}
Eficiência = {ribo*100:.1f}%"""
    
    ax.text(0.02, 0.98, texto_info, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', 
            alpha=0.9, edgecolor='#D1D5DB'))
    
    ax.set_xticks(np.arange(0, horas + 2, 2 if horas <= 24 else horas/12))
    
    buf = BytesIO()
    fig.savefig(buf, format='png', dpi=120, bbox_inches='tight', 
                facecolor='#F9FAFB', edgecolor='none')
    plt.close(fig)
    
    od_final = od[-1]
    if taxa > 0.01:
        tempo_duplicacao = 0.693 / taxa
    else:
        tempo_duplicacao = float('inf')
    
    return base64.b64encode(buf.getbuffer()).decode("ascii"), od_final, tempo_duplicacao, estado

@app.route('/', methods=['GET', 'POST'])
def index():
    opcoes_meio = list(simulador.meios.keys())
    meio = 'YPD (Yeast Extract Peptone Dextrose)'
    ph = 6.5
    temp = 37.0
    farnesol = 0.0
    od_inicial = 0.05
    horas = 48
    
    image_data = None
    diagnosticos = []
    scores = {}
    error_msg = None
    pop_final = 0.0
    t_dup = 0.0
    estado = ""
    
    if request.method == 'POST':
        try:
            meio = request.form.get('meio', meio)
            ph = float(request.form.get('ph', ph))
            temp = float(request.form.get('temperatura', temp))
            farnesol = float(request.form.get('farnesol', farnesol))
            od_inicial = float(request.form.get('od_inicial', od_inicial))
            horas = float(request.form.get('horas', horas))
            
            param = simulador.parametros
            
            if temp >= param['temperatura']['letal']:
                error_msg = "❌ MORTE: Temperatura letal (>55°C)"
            elif temp <= 0:
                error_msg = "⏸️ ESTASE: Temperatura de congelamento"
            elif ph <= param['ph']['letal_min'] or ph >= param['ph']['letal_max']:
                error_msg = "❌ MORTE: pH letal"
            elif temp < param['temperatura']['min']:
                error_msg = f"⚠️ Crescimento muito lento: Temperatura abaixo de {param['temperatura']['min']}°C"
            elif temp > param['temperatura']['max']:
                error_msg = f"⚠️ Estresse térmico severo: Temperatura acima de {param['temperatura']['max']}°C"
            elif ph < param['ph']['min']:
                error_msg = f"⚠️ Ambiente ácido extremo: pH abaixo de {param['ph']['min']}"
            elif ph > param['ph']['max']:
                error_msg = f"⚠️ Ambiente alcalino extremo: pH acima de {param['ph']['max']}"
            else:
                scores, od_max = simulador.prever(meio, ph, temp, farnesol)
                ribo = scores['Atividade Ribossomal']
                image_data, pop_final, t_dup, estado = gerar_grafico(ribo, od_inicial, horas, od_max)

                try:
                    print('\n=== Simulação executada ===')
                    print(f'Meio: {meio} | pH: {ph} | Temp: {temp}°C | Farnesol: {farnesol}')
                    print(f'OD inicial: {od_inicial} | Horas: {horas}')
                    print(f'OD_max (meio): {od_max}')
                    print(f'OD_final (simulação): {pop_final}')
                    print(f'Tempo de duplicação (t_dup): {t_dup}')
                    print(f'Atividade Ribossomal (ribo): {ribo}')
                    print('Scores:')
                    for k, v in scores.items():
                        try:
                            print(f'  {k}: {v:.6f}')
                        except Exception:
                            print(f'  {k}: {v}')
                    print('===========================\n')
                except Exception as e:
                    print('Erro ao imprimir resultados:', e)
                
                diagnosticos = []
                meio_info = simulador.meios.get(meio)
                
                if ribo > 0.8:
                    diagnosticos.append("✅ Crescimento ótimo - Fase exponencial ativa")
                elif ribo > 0.6:
                    diagnosticos.append("📈 Crescimento bom - Taxa elevada")
                elif ribo > 0.4:
                    diagnosticos.append("📊 Crescimento moderado - Condições sub-ótimas")
                elif ribo > 0.2:
                    diagnosticos.append("🐌 Crescimento lento - Fatores limitantes")
                elif ribo > 0.05:
                    diagnosticos.append("⚠️ Crescimento mínimo - Próximo à fase estacionária")
                else:
                    diagnosticos.append("⏸️ Latência/Declínio - Crescimento insignificante")
                
                hifas = scores['Formação de Hifas']
                if hifas > 0.6:
                    diagnosticos.append("🔄 Transição completa para hifas - Virulência elevada")
                elif hifas > 0.4:
                    diagnosticos.append("🌱 Formação de pseudohifas - Transição ativa")
                elif hifas > 0.2:
                    diagnosticos.append("🟡 Início de filamentação - Primeiros sinais")
                elif hifas > 0.05:
                    diagnosticos.append("🔵 Forma mista - Levedura predominante")
                else:
                    diagnosticos.append("⚪ Forma exclusiva de levedura")
                
                if farnesol > 5.0:
                    if hifas < 0.3:
                        diagnosticos.append("🧪 Farnesol ativo: inibição significativa de hifas")
                    else:
                        diagnosticos.append("⚠️ Farnesol presente mas com efeito limitado")
                
                hog1 = scores['Via MAPK Hog1']
                if hog1 > 0.7:
                    diagnosticos.append("🔥 Estresse severo: ativação máxima da via Hog1")
                elif hog1 > 0.5:
                    diagnosticos.append("⚠️ Estresse moderado: resposta compensatória ativa")
                elif hog1 > 0.3:
                    diagnosticos.append("📊 Estresse leve: mecanismos adaptativos ativados")
                elif hog1 > 0.1:
                    diagnosticos.append("🔬 Resposta basal: via Hog1 levemente ativada")
                else:
                    diagnosticos.append("🟢 Sem estresse: via Hog1 inativa")
                
                if temp > 39.0:
                    diagnosticos.append(f"🌡️ Hipertermia: {temp:.1f}°C (acima da ótima)")
                elif temp < 30.0:
                    diagnosticos.append(f"❄️ Hipotermia: {temp:.1f}°C (abaixo da ótima)")
                
                diff_ph = abs(ph - meio_info['ph_otimo']) if meio_info else abs(ph - 6.5)
                if diff_ph > 1.0:
                    diagnosticos.append(f"⚗️ pH desviado: {ph:.1f} vs ótimo {meio_info['ph_otimo']:.1f}")
                
                if meio_info:
                    diagnosticos.append(f"🧫 Meio: {meio_info['descricao']}")
                
                if t_dup < float('inf'):
                    if t_dup < 1.5:
                        diagnosticos.append(f"⚡ Duplicação rápida: {t_dup:.1f} horas")
                    elif t_dup < 3.0:
                        diagnosticos.append(f"📊 Duplicação moderada: {t_dup:.1f} horas")
                    else:
                        diagnosticos.append(f"🐌 Duplicação lenta: {t_dup:.1f} horas")
                
        except Exception as e:
            error_msg = f"❌ Erro técnico: {str(e)}"
    
    return render_template('index.html',
                           opcoes_meio=opcoes_meio,
                           meio=meio,
                           ph=ph,
                           temp=temp,
                           farnesol=farnesol,
                           od_inicial=od_inicial,
                           horas=horas,
                           image_data=image_data,
                           diagnosticos=diagnosticos,
                           scores=scores,
                           error_msg=error_msg,
                           pop_final=pop_final,
                           t_dup=t_dup,
                           estado=estado,
                           parametros=simulador.parametros,
                           simulador=simulador,
                           float=float)

if __name__ == '__main__':
    port = int(os.environ.get('PORT', 5000))
    print("\n" + "="*70)
    print("🔬 SIMULADOR  DE CANDIDA ALBICANS")
    print("="*70)
    print("📚 Baseado em literatura científica revisada:")
    print("   • Davis et al., 2000 - pH ótimo e faixas de crescimento")
    print("   • Hornby et al., 2001 - Efeito do farnesol")
    print("   • Lachke et al., 2002 - Cinética de crescimento")
    print("   • Lorenz et al., 2000 - Formação de hifas")
    print("="*70)
    print(f"🌐 Acesse: http://localhost:{port}")
    print("="*70)
    app.run(host='0.0.0.0', port=port, debug=True, use_reloader=False)