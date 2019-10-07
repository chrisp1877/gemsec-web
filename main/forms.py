from django import forms

class MD_Visualization_Form(forms.Form):
    analysis_choices = [("cos_sin_pca", "Cos/Sin PCA"), ("hilbert", "Hilbert PCA"), ("double_hilbert", "Double Hilbert PCA")]
    viz_choices = [("prob_map", "Probability Map"), ("density_plot", "Density Plot"), ("pca_gif", "PCA Gif")]
    analysis = forms.ChoiceField(choices=analysis_choices, widget=forms.RadioSelect)
    viz = forms.ChoiceField(choices=viz_choices, widget=forms.RadioSelect)
    file = forms.FileField()

    
'''
class MD_Visualization_Full_Form(forms.modelForm):
    class meta:
        pass
'''