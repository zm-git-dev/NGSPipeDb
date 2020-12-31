# Generated by Django 3.1.4 on 2020-12-31 05:28

from django.db import migrations, models


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Autoincrements',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('base', models.TextField(blank=True, null=True)),
                ('n', models.IntegerField(blank=True, null=True)),
            ],
            options={
                'db_table': 'autoincrements',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='Directives',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('directive', models.TextField(blank=True, null=True)),
            ],
            options={
                'db_table': 'directives',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='Duplicates',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('idspecid', models.TextField(blank=True, null=True)),
                ('newid', models.TextField(blank=True, null=True)),
            ],
            options={
                'db_table': 'duplicates',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='Features',
            fields=[
                ('id', models.TextField(blank=True, primary_key=True, serialize=False)),
                ('seqid', models.TextField(blank=True, null=True)),
                ('source', models.TextField(blank=True, null=True)),
                ('featuretype', models.TextField(blank=True, null=True)),
                ('start', models.IntegerField(blank=True, null=True)),
                ('end', models.IntegerField(blank=True, null=True)),
                ('score', models.TextField(blank=True, null=True)),
                ('strand', models.TextField(blank=True, null=True)),
                ('frame', models.TextField(blank=True, null=True)),
                ('attributes', models.TextField(blank=True, null=True)),
                ('extra', models.TextField(blank=True, null=True)),
                ('bin', models.IntegerField(blank=True, null=True)),
            ],
            options={
                'db_table': 'features',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='Meta',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('dialect', models.TextField(blank=True, null=True)),
                ('version', models.TextField(blank=True, null=True)),
            ],
            options={
                'db_table': 'meta',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='Relations',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('parent', models.TextField(blank=True, null=True)),
                ('child', models.TextField(blank=True, null=True)),
                ('level', models.IntegerField(blank=True, null=True)),
            ],
            options={
                'db_table': 'relations',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='SqliteStat1',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('tbl', models.TextField(blank=True, null=True)),
                ('idx', models.TextField(blank=True, null=True)),
                ('stat', models.TextField(blank=True, null=True)),
            ],
            options={
                'db_table': 'sqlite_stat1',
                'managed': False,
            },
        ),
    ]
