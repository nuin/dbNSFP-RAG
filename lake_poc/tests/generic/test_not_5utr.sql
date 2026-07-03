{% test not_5utr(model, column_name) %}
select {{ column_name }}
from {{ model }}
where regexp_matches({{ column_name }}, '^c\.-\d')
{% endtest %}
