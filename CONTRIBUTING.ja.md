# コントリビューションガイド

> English: [CONTRIBUTING.md](CONTRIBUTING.md)(原本)

このファイルは日本語訳です。原本は `CONTRIBUTING.md` です。変更するときは、同じコミットで両方を更新してください。

## ブランチモデル(GitFlow)

このプロジェクトは GitFlow 形式のブランチモデルに従います。

- **`main`**: 最新のリリース済みの状態です。ここへ直接コミットしません。`release/*` と `hotfix/*` からのマージだけを受け付け、マージごとにバージョン番号でタグを付けます(例: `2026.10.0`。`v` は付けません)。
- **`develop`**: 開発を統合するブランチです。日々の作業はこのブランチを対象にします。
- **`feature/*`、`fix/*`、`chore/*`**: 作業単位ごとに `develop` から分岐し、レビュー後に `develop` へマージします。
- **`release/X.Y.Z`**: リリースを準備するときに `develop` から分岐します。版の更新などの最終調整をここで行い、`main`(タグ付き)と `develop` の両方にマージします。
- **`hotfix/*`**: リリース済みの版への緊急の修正のため、`main` から分岐します。`main`(タグ付き)と `develop` の両方にマージします。

作業ブランチごとに git worktree を持ち、`develop` の worktree の隣に作成されます(`devtool/flow.sh start <branch>`)。

## バージョン

バージョンは `YYYY.M.PATCH`(暦に基づく版番号)に従います。ある年月の最初のリリースは `PATCH=0`(例: `2026.10.0`)です。同じ月の後続のリリースでは `PATCH` を増やします。

バージョンは `qclobot/__version__.py` の `__version_info__` タプルに保存され、`setup.cfg` がこれを読みます。`devtool/flow.sh start release/<version>` と `hotfix/<version>` がこれを更新します。

## テストと CI の実行

継続的インテグレーションは `scripts/ci-local.sh` でローカルに実行します。

### 前提

- Python 3.9 以降で、`venv` モジュールが使えること(`PYTHON` で使うインタプリタを指定できます。既定は `python3`)。
- そのインタプリタに ProteinDF_bridge がインストールされていること。テストは `proteindf_bridge` を import します。
- AmberTools がインストールされ、`AMBERHOME` がそれを指していること。`unit` ジョブは `$AMBERHOME/bin/tleap`、`sander`、`ambpdb`、`reduce` を必要とします。テストは `AMBERHOME` を通してこれらを呼び出します。

### ローカル CI の実行

```bash
scripts/ci-local.sh              # 既定のジョブをすべて実行: unit, docs
scripts/ci-local.sh docs         # 1つのジョブだけ実行
scripts/ci-local.sh --help       # 使い方
```

ジョブ:

- `unit`: 作業ツリーの複製を一時ディレクトリに作り、その上で `tests/test_*.py` を `unittest` で実行します。テストの出力は worktree に現れません。`AMBERHOME` または AmberTools の実行ファイルが無い場合、または `proteindf_bridge` を import できない場合は、テストの前に失敗します。前提が無いことはスキップではなく失敗として扱います。
- `docs`: 作業ツリーの複製から、英語と日本語の Sphinx HTML を隔離した仮想環境でビルドします。

## ドキュメント(英語と日本語)

英語が原本です。ドキュメントを変更するときは、日本語訳を同じコミットに含めてください。

- Sphinx のソース(`doc/source/*.rst`): 英語のソースを変更したら、`doc/source/locale/ja/LC_MESSAGES/` の日本語カタログを更新します(`.pot` を再生成して `sphinx-intl update` を実行し、新しいメッセージを翻訳します)。
- Markdown のドキュメント(`*.md`): 英語のファイルの隣に `*.ja.md` を置き、両方を更新します。

`docs` ジョブは両方の言語をビルドし、どちらかが失敗すると失敗します。

## API リファレンス(`doc/source`)

`qclobot/` 配下のモジュールを追加・削除・改名したら、`doc/source/` の Sphinx API ページを更新し、`docs` ジョブで文書をビルドし直してください。

## AI エージェントとの開発

| 担当 | 役割 | 規則 |
|---|---|---|
| ユーザー | 方向性と設計判断、マージの承認 | |
| Claude | タスク指示書の作成、レビュー、承認後のマージ | `CLAUDE.md` |
| agy / GitHub Copilot | 作業ブランチ専用の worktree での実装 | `AGENTS.md` |

1. **タスク指示書**: Claude が `TASK_<name>.md` を作成します。作業ブランチ(`**Branch**:` の行)、対象、完了の定義を書きます。タスク指示書は非公開です。非公開リポジトリ `ProteinDF_private` の `qclobot/tasks/` にあり、`develop` の worktree から `private/` として参照します(git では無視されます)。
2. **実装**: タスクをエージェントに渡します。ブランチと worktree が無ければ、`develop` の worktree の隣に作成されます。
   ```bash
   devtool/delegate.sh private/tasks/TASK_<name>.md            # agy(非対話、自動承認)
   devtool/delegate.sh -a copilot private/tasks/TASK_<name>.md # GitHub Copilot CLI
   devtool/delegate.sh -i private/tasks/TASK_<name>.md         # 対話モード
   ```
   出力は `.git/agent-logs/` に保存されます。エージェントがアクセスしてよい追加のディレクトリ(ローカルの ProteinDF_bridge のチェックアウトや AmberTools の場所など)は、`$(git rev-parse --git-common-dir)/delegate.conf` に1行ずつ書きます。エージェントが利用上限で止まると、スクリプトは終了コード 75 で終わります。`-c` で再開するか、`-W` で待ちます。
3. **レビュー**: Claude が `devtool/review.sh <branch>` を実行します。これは `develop` 側の `scripts/ci-local.sh` を作業ブランチの worktree に対して実行し、結果をタスク指示書に記録します。
4. **マージ**: ユーザーの承認後に `devtool/flow.sh finish <branch>` を実行します。リリースは `devtool/flow.sh start release/<YYYY.M.PATCH>`(版を更新)と `devtool/flow.sh finish release/<YYYY.M.PATCH>`(`main` へのマージ、タグ付け、`develop` への取り込み)を使います。

`devtool/tasks.sh` は、マージされていないタスクと作業ブランチの一覧を表示します(`--all` でマージ済みも含みます)。これらのスクリプトは何も push しません。

## 公開と非公開の情報

公開するファイル(追跡されるファイル、コミットメッセージ、ドキュメント)に、ローカルのディレクトリのパス、タスク指示書の名前、非公開リポジトリの内容を書いてはいけません。ローカルの場所は環境変数か、`$(git rev-parse --git-common-dir)/delegate.conf` に置きます。
