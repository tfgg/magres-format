{% extends 'markdown.tpl' %}

{% block input %}
```python
{{ cell.input}}
```
{% endblock input %}


{% block stream %}
<div class='stream'>
<pre>{{ output.text|escape }}</pre>
</div>
{% endblock %}


{% block pyout %}
<div class='stream'>
<pre>{{ output.text|escape }}</pre>
</div>
{% endblock %}
