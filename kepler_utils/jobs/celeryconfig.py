from kombu import Exchange, Queue

BROKER_URL = 'amqp://guest@localhost//'
CELERY_RESULT_BACKEND = 'amqp://'

CELERY_TASK_SERIALIZER = 'json'
CELERY_RESULT_SERIALIZER = 'json'
CELERY_ACCEPT_CONTENT=['json']
CELERY_TIMEZONE = 'US/Pacific'
CELERY_ENABLE_UTC = True

CELERYD_CONCURRENCY = 1

CELERYD_PREFETCH_MULTIPLIER = 1

CELERY_DEFAULT_QUEUE = 'default'
CELERY_QUEUES = (
    Queue ('default', Exchange ('default'), routing_key = 'default'),
    Queue ('priority', Exchange ('priority'), routing_key = 'priority'),
)